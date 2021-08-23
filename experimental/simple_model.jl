
function w_partial_t(w_x, w_y, phi_x, phi_y, psi_x, psi_y, am)

    _a = -(phi_x.*w_y - phi_y.*w_x)
    _b = -(psi_x.*w_y - psi_y.*w_x)
    return _a + am*_b
end


function psi_partial_t(psi_x, psi_y, phi_x, phi_y)

    return -(phi_x.*psi_y - phi_y.*psi_x)
end


function forward_euler(x, x_t, K1, K2, dt, gd::GhostData)

    return (x + dt*(K1*x_t)) ./ (1 .+ dt*diag(K2))
end


function solve_pde(A, x0, b, gd::GhostData)

    An, _ = require_boundary_conditions(A, gd)
    bn = gd.Proj*b

    bn = transpose(An)*bn
    An = transpose(An)*An
    # x = jacobi_preconditioned_gmres(An, gd.Proj*x0, bn, 1, 1)
    x = An \ bn
    return gd.R*transpose(gd.Proj)*x
end


function solve_diffusion(x, k, Dxx, Dyy, gd::GhostData)

    M = I - k*(Dxx + Dyy)
    return solve_pde(M, x, x, gd)
end


function solve_vorticity_eqn(phi, w, Dxx, Dyy, gd::GhostData)

    return solve_pde(Dxx + Dyy, phi, w, gd)
end


function leapfrog!(w, psi, phi, Dx, Dxx, Dy, Dyy, K1, K2, gd, dt, am, k_psi, k_w)

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
    psi[:,3] = forward_euler(0.5*(psi[:,1] + psi[:,2]), psi_t, K1, K2, dt, gd)
    w[:,3] = forward_euler(0.5*(w[:,1] + w[:,2]), w_t, K1, K2, dt, gd)

    # compute diffusion
    psi[:,3] = solve_diffusion(psi[:,3], k_psi, Dxx, Dyy, gd)
    w[:,3] = solve_diffusion(w[:,3], k_w, Dxx, Dyy, gd)

    # solve vorticity equation
    phi[:] = solve_vorticity_eqn(phi, w[:,3], Dxx, Dyy, gd)

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
    psi[:,3] = forward_euler(psi[:,2], psi_t, K1, K2, dt, gd)
    w[:,3] = forward_euler(w[:,2], w_t, K1, K2, dt, gd)

    # compute diffusion
    psi[:,3] = solve_diffusion(psi[:,3], k_psi, Dxx, Dyy, gd)
    w[:,3] = solve_diffusion(w[:,3], k_w, Dxx, Dyy, gd)

    # solve vorticity equation
    phi[:] = solve_vorticity_eqn(phi, w[:,3], Dxx, Dyy, gd)

    # shift indices
    psi[:,1:2] = psi[:,2:3]
    w[:,1:2] = w[:,2:3]
end
