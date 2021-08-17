
function w_partial_t(w_x, w_y, phi_x, phi_y, psi_x, psi_y, am)

    _a = -(phi_x.*w_y - phi_y.*w_x)
    _b = -(psi_x.*w_y - psi_y.*w_x)
    return _a + am*_b
end


function psi_partial_t(psi_x, psi_y, phi_x, phi_y)

    return -(phi_x.*psi_y - phi_y.*psi_x)
end


function forward_euler(x, x_t, K1, K2, dt, gc::GhostConditions)

    return (x + dt*(K1*x_t)) ./ (1 .+ dt*diag(K2))
end


function solve_pde(A, x0, b, gc::GhostConditions)

    An, _ = require_boundary_conditions(A, gc)
    bn = gc.Proj*b
    # x = jacobi_preconditioned_gmres(An, gc.Proj*x0, bn, 1, 1)
    x = An \ bn
    return gc.Mirror*transpose(gc.Proj)*x
end


function solve_diffusion(x, k, Dxx, Dyy, gc::GhostConditions)

    M = I - k*(Dxx + Dyy)
    return solve_pde(M, x, x, gc)
end


function solve_vorticity_eqn(phi, w, Dxx, Dyy, gc::GhostConditions)

    return solve_pde(Dxx + Dyy, phi, w, gc)
end


function leapfrog!(w, psi, phi, Dx, Dxx, Dy, Dyy, K1, K2, gc, dt, am, k_psi, k_w)

    # take derivatives
    psi_x = Dx*psi[:,2]
    psi_y = Dy*psi[:,2]
    phi_x = Dx*phi
    phi_y = Dy*phi
    w_x = Dx*w[:,2]
    w_y = Dy*w[:,2]

    w_t = w_partial_t(w_x, w_y, phi_x, phi_y, psi_x, psi_y, am)
    psi_t = psi_partial_t(psi_x, psi_y, phi_x, phi_y)

    # compute predictor step
    psi[:,3] = forward_euler(0.5*(psi[:,1] + psi[:,2]), psi_t, K1, K2, dt, gc)
    w[:,3] = forward_euler(0.5*(w[:,1] + w[:,2]), w_t, K1, K2, dt, gc)

    # compute diffusion
    psi[:,3] = solve_diffusion(psi[:,3], k_psi, Dxx, Dyy, gc)
    w[:,3] = solve_diffusion(w[:,3], k_w, Dxx, Dyy, gc)

    # solve vorticity equation
    phi[:] = solve_vorticity_eqn(phi, w[:,3], Dxx, Dyy, gc)

    # take derivatives again
    psi_x = Dx*psi[:,3]
    psi_y = Dy*psi[:,3]
    phi_x = Dx*phi
    phi_y = Dy*phi
    w_x = Dx*w[:,3]
    w_y = Dy*w[:,3]

    w_t = w_partial_t(w_x, w_y, phi_x, phi_y, psi_x, psi_y, am)
    psi_t = psi_partial_t(psi_x, psi_y, phi_x, phi_y)

    # compute corrector step
    psi[:,3] = forward_euler(psi[:,2], psi_t, K1, K2, dt, gc)
    w[:,3] = forward_euler(w[:,2], w_t, K1, K2, dt, gc)

    # compute diffusion
    psi[:,3] = solve_diffusion(psi[:,3], k_psi, Dxx, Dyy, gc)
    w[:,3] = solve_diffusion(w[:,3], k_w, Dxx, Dyy, gc)

    # solve vorticity equation
    phi[:] = solve_vorticity_eqn(phi, w[:,3], Dxx, Dyy, gc)

    # shift indices
    psi[:,1:2] = psi[:,2:3]
    w[:,1:2] = w[:,2:3]
end
