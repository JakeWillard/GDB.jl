

function bval_u(Te, Ti, trgt_sgn, H3)

    cs = trgt_sgn .* sqrt.(Te + Ti)

    return H3 * cs
end


function bval_phi(phi, Te, lcfs_avg, H1, H2, H3)

    lambda = 2.695
    phi_in = dot(lcfs_avg, phi)*ones(length(phi))
    phi_out = 2.695 * Te

    return H1*phi_in + (H2 + H3)*phi_out
end


function solve_pde(A, x0, b, xb, gc::GhostConditions)

    An, B = require_boundary_conditions(A, gc)
    bn = gc.Proj*b - B*xb
    x = jacobi_preconditioned_gmres(An, gc.Proj*x0, bn, 1, 1)
    # x = An \ bn
    return gc.Mirror*transpose(gc.Proj)*x + (I - gc.Mirror)*xb
end


function solve_diffusion(x, xb, k, Dxx, Dyy, gc::GhostConditions)

    M = I - k*(Dxx + Dyy)
    return solve_pde(M, x, x, xb, gc)
end


function solve_vorticity_eqn(phi, phi_b, w, n, lnn_x, lnn_y, Pi_xx, Pi_yy, ad, Dx, Dy, Dxx, Dyy, gc::GhostConditions)

    M = Dxx + Dyy + Diagonal(lnn_x) * Dx + Diagonal(lnn_y) * Dy
    b = (w - ad*(Pi_xx + Pi_yy)) ./ n

    return solve_pde(M, phi, b, phi_b, gc)
end


function solve_helmholtz_eqn(psi, A, de2, Dxx, Dyy, gc::GhostConditions)

    M = I - de2*(Dxx + Dyy)
    return solve_pde(M, psi, A, zeros(length(A)), gc)
end


function forward_euler(x, x_t, xb, K1, K2, dt, gc::GhostConditions)

    return (x + dt*(K1*x_t + K2*gc.Mirror*xb)) ./ (1 .+ dt*diag(K2))
end


function leapfrog!(lnn, lnTe, lnTi, u, w, A, phi, psi, n, Te, Ti, Pe, Pi, j, jn, Sn, STe, STi, Dx, Dy, Dxy, Dxx, Dyy, Dxxx, Dyyy,
                   Dxxy, Dxyy, Ds, Dss, dt, N_subcycle, K1, K2, H1, H2, H3, trgt_sgn, lcfs_avg, GC_nmann, GC_dchlt, GC_u, kdiff_lnn,
                   kdiff_lnTe, kdiff_lnTi, kdiff_u, kdiff_w, kdiff_A, am, ad, ki, ke, er, eg, ev, de2, eta)

    # homogeneous boundary values
    bval_hom = zeros(length(Te))

    # take derivatives of psi
    psi_x = Dx*psi
    psi_y = Dy*psi
    psi_xy = Dxy*psi
    psi_xx = Dxx*psi
    psi_yy = Dyy*psi

    # subcycle electron thermal conduction terms
    # lnTe0 = lnTe[:,2]
    # for _=1:N_subcycle
    #
    #     # derivatives
    #     lnTe_x = Dx*lnTe[:,2]
    #     lnTe_y = Dy*lnTe[:,2]
    #     lnTe_xy = Dxy*lnTe[:,2]
    #     lnTe_xx = Dxx*lnTe[:,2]
    #     lnTe_yy = Dyy*lnTe[:,2]
    #     lnTe_s = partial_s(lnTe[:,2], lnTe_x, lnTe_y, psi_x, psi_y, Ds, am)
    #     lnTe_ss = partial_ss(lnTe[:,2], lnTe_s, lnTe_x, lnTe_y, lnTe_xy, lnTe_xx, lnTe_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, Dss, Ds, Dx, Dy, am)
    #
    #     lnTe_t = lnT_partial_thermal(lnTe_s, lnTe_ss, Te, n, ke) + STe
    #     lnTe[:,3] = forward_euler(0.5*(lnTe[:,1] + lnTe[:,2]), lnTe_t, bval_hom, K1, K2, dt/N_subcycle, GC_nmann)
    #
    #     # derivatives again
    #     lnTe_x = Dx*lnTe[:,2]
    #     lnTe_y = Dy*lnTe[:,2]
    #     lnTe_xy = Dxy*lnTe[:,2]
    #     lnTe_xx = Dxx*lnTe[:,2]
    #     lnTe_yy = Dyy*lnTe[:,2]
    #     lnTe_s = partial_s(lnTe[:,2], lnTe_x, lnTe_y, psi_x, psi_y, Ds, am)
    #     lnTe_ss = partial_ss(lnTe[:,2], lnTe_s, lnTe_x, lnTe_y, lnTe_xy, lnTe_xx, lnTe_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, Dss, Ds, Dx, Dy, am)
    #
    #     lnTe_t = lnT_partial_thermal(lnTe_s, lnTe_ss, Te, n, ke) + STe
    #     lnTe[:,3] = forward_euler(lnTe[:,2], lnTe_t, bval_hom, K1, K2, dt/N_subcycle, GC_nmann)
    #
    #     # shift indices
    #     lnTe[:,1:2] = lnTe[:,2:3]
    # end
    # lnTe[:,2] = 0.5*(lnTe0 + lnTe[:,2])
    #
    # # subcycle ion thermal conduction terms
    # lnTi0 = lnTi[:,2]
    # for _=1:N_subcycle
    #
    #     # derivatives
    #     lnTi_x = Dx*lnTi[:,2]
    #     lnTi_y = Dy*lnTi[:,2]
    #     lnTi_xy = Dxy*lnTi[:,2]
    #     lnTi_xx = Dxx*lnTi[:,2]
    #     lnTi_yy = Dyy*lnTi[:,2]
    #     lnTi_s = partial_s(lnTi[:,2], lnTi_x, lnTi_y, psi_x, psi_y, Ds, am)
    #     lnTi_ss = partial_ss(lnTi[:,2], lnTi_s, lnTi_x, lnTi_y, lnTi_xy, lnTi_xx, lnTi_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, Dss, Ds, Dx, Dy, am)
    #
    #     lnTi_t = lnT_partial_thermal(lnTi_s, lnTi_ss, Te, n, ke) + STi
    #     lnTi[:,3] = forward_euler(0.5*(lnTi[:,1] + lnTi[:,2]), lnTi_t, bval_hom, K1, K2, dt/N_subcycle, GC_nmann)
    #
    #     # derivatives again
    #     lnTi_x = Dx*lnTi[:,2]
    #     lnTi_y = Dy*lnTi[:,2]
    #     lnTi_xy = Dxy*lnTi[:,2]
    #     lnTi_xx = Dxx*lnTi[:,2]
    #     lnTi_yy = Dyy*lnTi[:,2]
    #     lnTi_s = partial_s(lnTi[:,2], lnTi_x, lnTi_y, psi_x, psi_y, Ds, am)
    #     lnTi_ss = partial_ss(lnTi[:,2], lnTi_s, lnTi_x, lnTi_y, lnTi_xy, lnTi_xx, lnTi_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, Dss, Ds, Dx, Dy, am)
    #
    #     lnTi_t = lnT_partial_thermal(lnTi_s, lnTi_ss, Te, n, ke) + STi
    #     lnTi[:,3] = forward_euler(lnTi[:,2], lnTi_t, bval_hom, K1, K2, dt/N_subcycle, GC_nmann)
    #
    #     # shift indices
    #     lnTi[:,1:2] = lnTi[:,2:3]
    # end
    # lnTi[:,2] = 0.5*(lnTi0 + lnTi[:,2])

    # take other derivatives
    phi_x = Dx*phi
    phi_y = Dy*phi
    phi_xy = Dxy*phi
    phi_xx = Dxx*phi
    phi_yy = Dyy*phi
    phi_xxx = Dxxx*phi
    phi_yyy = Dyyy*phi
    phi_xxy = Dxxy*phi
    phi_xyy = Dxyy*phi
    phi_c = -2*phi_y[:]
    phi_s = partial_s(phi, phi_x, phi_y, psi_x, psi_y, Ds, am)

    Pe_x = Dx*Pe
    Pe_y = Dy*Pe
    Pe_xy = Dxy*Pe
    Pe_xx = Dxx*Pe
    Pe_yy = Dyy*Pe
    Pe_xxx = Dxxx*Pe
    Pe_yyy = Dyyy*Pe
    Pe_xxy = Dxxy*Pe
    Pe_xyy = Dxyy*Pe
    Pe_c = -2*Pe_y[:]
    Pe_s = partial_s(Pe, Pe_x, Pe_y, psi_x, psi_y, Ds, am)

    Pi_x = Dx*Pi
    Pi_y = Dy*Pi
    Pi_c = -2*Pi_y[:]
    Pi_s = partial_s(Pi, Pi_x, Pi_y, psi_x, psi_y, Ds, am)

    lnn_x = Dx*lnn[:,2]
    lnn_y = Dy*lnn[:,2]
    lnn_xy = Dxy*lnn[:,2]
    lnn_xx = Dxx*lnn[:,2]
    lnn_yy = Dyy*lnn[:,2]
    lnn_s = partial_s(lnn[:,2], lnn_x, lnn_y, psi_x, psi_y, Ds, am)

    lnTe_x = Dx*lnTe[:,2]
    lnTe_y = Dy*lnTe[:,2]
    lnTe_c = -2*lnTe_y[:]
    lnTe_s = partial_s(lnTe[:,2], lnTe_x, lnTe_y, psi_x, psi_y, Ds, am)

    lnTi_x = Dx*lnTi[:,2]
    lnTi_y = Dy*lnTi[:,2]
    lnTi_c = -2*lnTi_y[:]
    lnTi_s = partial_s(lnTi[:,2], lnTi_x, lnTi_y, psi_x, psi_y, Ds, am)

    u_x = Dx*u[:,2]
    u_y = Dy*u[:,2]
    u_c = -2*u_y[:]
    u_s = partial_s(u[:,2], u_x, u_y, psi_x, psi_y, Ds, am)

    j_x = Dx*j
    j_y = Dy*j
    j_s = partial_s(j, j_x, j_y, psi_x, psi_y, Ds, am)

    jn_x = Dx*jn
    jn_y = Dy*jn
    jn_s = partial_s(jn, jn_x, jn_y, psi_x, psi_y, Ds, am)

    G = G_def(Ti, phi_c, Pi_c, u_s, er, ev, ad)
    G_x = Dx*G
    G_y = Dy*G
    G_c = -2*G_y[:]
    G_s = partial_s(G, G_x, G_y, psi_x, psi_y, Ds, am)

    # evaluate total derivatives
    lnn_Dt = lnn_total_ion_derivative(phi_c, Pe_c, n, j_s, u_s, er, ev, ad)
    lnTe_Dt = lnTe_total_electron_derivative_adiabatic(lnn_Dt, j, lnn_s, n, Te, lnTe_c, er, ad)
    lnTi_Dt = lnTi_total_ion_derivative_adiabatic(lnn_Dt, Ti, lnTi_c, er, ad)
    u_Dt = u_total_ion_derivative(Pe_s, Pi_s, G_s, n, Ti, u_c, ev, eg, er, ad)

    # evaluate partial derivatives
    lnn_t = general_partial_derivative(lnn_Dt, lnn_x, lnn_y, lnn_s, phi_x, phi_y, u[:,2]) + Sn
    lnTe_t = general_partial_derivative(lnTe_Dt, lnTe_x, lnTe_y, lnTe_s, phi_x, phi_y, u[:,2] - jn / ev) + STe
    lnTi_t = general_partial_derivative(lnTi_Dt, lnTi_x, lnTi_y, lnTi_s, phi_x, phi_y, u[:,2]) + STi
    u_t = general_partial_derivative(u_Dt, u_x, u_y, u_s, phi_x, phi_y, u[:,2])
    w_t = w_partial_t(n, lnn_x, lnn_y, lnn_xy, lnn_xx, lnn_yy,
                         Pe_x, Pe_y, Pe_xy, Pe_xx, Pe_yy, Pe_xxx, Pe_yyy, Pe_xxy, Pe_xyy, Pe_c,
                         phi_x, phi_y, phi_xy, phi_xx, phi_yy, phi_xxx, phi_yyy, phi_xxy, phi_xyy,
                         G_c, Pi_c, j_s, ad, eg)
    A_t = A_partial_t(phi_x, phi_y, phi_s, j, jn, jn_x, jn_y, jn_s, Pe_s, u[:,2], n, Te, de2, ev, ad, am, eta)

    # compute inhomogeneous boundary values
    u_b = bval_u(Te, Ti, trgt_sgn, H3)
    phi_b = bval_phi(phi, Te, lcfs_avg, H1, H2, H3)

    # compute predictor step
    lnn[:,3] = forward_euler(0.5*(lnn[:,1] + lnn[:,2]), lnn_t, bval_hom, K1, K2, dt, GC_nmann)
    lnTe[:,3] = forward_euler(0.5*(lnTe[:,1] + lnTi[:,2]), lnTe_t, bval_hom, K1, K2, dt, GC_nmann)
    lnTi[:,3] = forward_euler(0.5*(lnTi[:,1] + lnTi[:,2]), lnTi_t, bval_hom, K1, K2, dt, GC_nmann)
    u[:,3] = forward_euler(0.5*(u[:,1] + u[:,2]), u_t, u_b, K1, K2, dt, GC_u)
    w[:,3] = forward_euler(0.5*(w[:,1] + w[:,2]), w_t, bval_hom, K1, K2, dt, GC_dchlt)
    A[:,3] = forward_euler(0.5*(A[:,1] + A[:,2]), A_t, bval_hom, K1, K2, dt, GC_dchlt)

    # solve diffusion problems
    lnn[:,3] = solve_diffusion(lnn[:,3], bval_hom, kdiff_lnn, Dxx, Dyy, GC_nmann)
    lnTe[:,3] = solve_diffusion(lnTe[:,3], bval_hom, kdiff_lnTe, Dxx, Dyy, GC_nmann)
    lnTi[:,3] = solve_diffusion(lnTi[:,3], bval_hom, kdiff_lnTi, Dxx, Dyy, GC_nmann)
    u[:,3] = solve_diffusion(u[:,3], u_b, kdiff_u, Dxx, Dyy, GC_u)
    w[:,3] = solve_diffusion(w[:,3], bval_hom, kdiff_w, Dxx, Dyy, GC_nmann)
    A[:,3] = solve_diffusion(A[:,3], bval_hom, kdiff_A, Dxx, Dyy, GC_nmann)

    # exponentiate log terms
    n[:] = exp.(lnn[:,3])
    Te[:] = exp.(lnTe[:,3])
    Ti[:] = exp.(lnTi[:,3])
    Pe[:] = n .* Te
    Pi[:] = n .* Ti

    # solve vorticity and helmholtz equations
    lnn_x = Dx * lnn[:,3]
    lnn_y = Dy * lnn[:,3]
    Pi_xx = Dxx * Pi
    Pi_yy = Dyy * Pi
    phi = solve_vorticity_eqn(phi, phi_b, w[:,3], n, lnn_x, lnn_y, Pi_xx, Pi_yy, ad, Dx, Dy, Dxx, Dyy, GC_dchlt)
    psi = solve_helmholtz_eqn(psi, A[:,3], de2, Dxx, Dyy, GC_dchlt)

    # compute current
    j[:] = (Dxx + Dyy) * psi
    jn[:] = j ./ n

    # take derivatives
    psi_x = Dx*psi
    psi_y = Dy*psi
    psi_xy = Dxy*psi
    psi_xx = Dxx*psi
    psi_yy = Dyy*psi

    phi_x = Dx*phi
    phi_y = Dy*phi
    phi_xy = Dxy*phi
    phi_xx = Dxx*phi
    phi_yy = Dyy*phi
    phi_xxx = Dxxx*phi
    phi_yyy = Dyyy*phi
    phi_xxy = Dxxy*phi
    phi_xyy = Dxyy*phi
    phi_c = -2*phi_y[:]
    phi_s = partial_s(phi, phi_x, phi_y, psi_x, psi_y, Ds, am)

    Pe_x = Dx*Pe
    Pe_y = Dy*Pe
    Pe_xy = Dxy*Pe
    Pe_xx = Dxx*Pe
    Pe_yy = Dyy*Pe
    Pe_xxx = Dxxx*Pe
    Pe_yyy = Dyyy*Pe
    Pe_xxy = Dxxy*Pe
    Pe_xyy = Dxyy*Pe
    Pe_c = -2*Pe_y[:]
    Pe_s = partial_s(Pe, Pe_x, Pe_y, psi_x, psi_y, Ds, am)

    Pi_x = Dx*Pi
    Pi_y = Dy*Pi
    Pi_c = -2*Pi_y[:]
    Pi_s = partial_s(Pi, Pi_x, Pi_y, psi_x, psi_y, Ds, am)

    lnn_x = Dx*lnn[:,3]
    lnn_y = Dy*lnn[:,3]
    lnn_xy = Dxy*lnn[:,3]
    lnn_xx = Dxx*lnn[:,3]
    lnn_yy = Dyy*lnn[:,3]
    lnn_s = partial_s(lnn[:,3], lnn_x, lnn_y, psi_x, psi_y, Ds, am)

    lnTe_x = Dx*lnTe[:,3]
    lnTe_y = Dy*lnTe[:,3]
    lnTe_c = -2*lnTe_y[:]
    lnTe_s = partial_s(lnTe[:,3], lnTe_x, lnTe_y, psi_x, psi_y, Ds, am)

    lnTi_x = Dx*lnTi[:,3]
    lnTi_y = Dy*lnTi[:,3]
    lnTi_c = -2*lnTi_y[:]
    lnTi_s = partial_s(lnTi[:,3], lnTi_x, lnTi_y, psi_x, psi_y, Ds, am)

    u_x = Dx*u[:,3]
    u_y = Dy*u[:,3]
    u_c = -2*u_y[:]
    u_s = partial_s(u[:,3], u_x, u_y, psi_x, psi_y, Ds, am)

    j_x = Dx*j
    j_y = Dy*j
    j_s = partial_s(j, j_x, j_y, psi_x, psi_y, Ds, am)

    jn_x = Dx*jn
    jn_y = Dy*jn
    jn_s = partial_s(jn, jn_x, jn_y, psi_x, psi_y, Ds, am)

    G = G_def(Ti, phi_c, Pi_c, u_s, er, ev, ad)
    G_x = Dx*G
    G_y = Dy*G
    G_c = -2*G_y[:]
    G_s = partial_s(G, G_x, G_y, psi_x, psi_y, Ds, am)

    # evaluate total derivatives
    lnn_Dt = lnn_total_ion_derivative(phi_c, Pe_c, n, j_s, u_s, er, ev, ad)
    lnTe_Dt = lnTe_total_electron_derivative_adiabatic(lnn_Dt, j, lnn_s, n, Te, lnTe_c, er, ad)
    lnTi_Dt = lnTi_total_ion_derivative_adiabatic(lnn_Dt, Ti, lnTi_c, er, ad)
    u_Dt = u_total_ion_derivative(Pe_s, Pi_s, G_s, n, Ti, u_c, ev, eg, er, ad)

    # evaluate partial derivatives
    lnn_t = general_partial_derivative(lnn_Dt, lnn_x, lnn_y, lnn_s, phi_x, phi_y, u[:,3])
    lnTe_t = general_partial_derivative(lnTe_Dt, lnTe_x, lnTe_y, lnTe_s, phi_x, phi_y, u[:,3] - jn / ev)
    lnTi_t = general_partial_derivative(lnTi_Dt, lnTi_x, lnTi_y, lnTi_s, phi_x, phi_y, u[:,3])
    u_t = general_partial_derivative(u_Dt, u_x, u_y, u_s, phi_x, phi_y, u[:,3])
    w_t = w_partial_t(n, lnn_x, lnn_y, lnn_xy, lnn_xx, lnn_yy,
                         Pe_x, Pe_y, Pe_xy, Pe_xx, Pe_yy, Pe_xxx, Pe_yyy, Pe_xxy, Pe_xyy, Pe_c,
                         phi_x, phi_y, phi_xy, phi_xx, phi_yy, phi_xxx, phi_yyy, phi_xxy, phi_xyy,
                         G_c, Pi_c, j_s, ad, eg)
    A_t = A_partial_t(phi_x, phi_y, phi_s, j, jn, jn_x, jn_y, jn_s, Pe_s, u[:,3], n, Te, de2, ev, ad, am, eta)

    # compute inhomogeneous boundary values
    u_b = bval_u(Te, Ti, trgt_sgn, H3)
    phi_b = bval_phi(phi, Te, lcfs_avg, H1, H2, H3)

    # compute corrector step
    lnn[:,3] = forward_euler(lnn[:,2], lnn_t, bval_hom, K1, K2, dt, GC_nmann)
    lnTe[:,3] = forward_euler(lnTi[:,2], lnTe_t, bval_hom, K1, K2, dt, GC_nmann)
    lnTi[:,3] = forward_euler(lnTi[:,2], lnTi_t, bval_hom, K1, K2, dt, GC_nmann)
    u[:,3] = forward_euler(u[:,2], u_t, u_b, K1, K2, dt, GC_u)
    w[:,3] = forward_euler(w[:,2], w_t, bval_hom, K1, K2, dt, GC_dchlt)
    A[:,3] = forward_euler(A[:,2], A_t, bval_hom, K1, K2, dt, GC_dchlt)

    # solve diffusion problems
    lnn[:,3] = solve_diffusion(lnn[:,3], bval_hom, kdiff_lnn, Dxx, Dyy, GC_nmann)
    lnTe[:,3] = solve_diffusion(lnTe[:,3], bval_hom, kdiff_lnTe, Dxx, Dyy, GC_nmann)
    lnTi[:,3] = solve_diffusion(lnTi[:,3], bval_hom, kdiff_lnTi, Dxx, Dyy, GC_nmann)
    u[:,3] = solve_diffusion(u[:,3], u_b, kdiff_u, Dxx, Dyy, GC_u)
    w[:,3] = solve_diffusion(w[:,3], bval_hom, kdiff_w, Dxx, Dyy, GC_nmann)
    A[:,3] = solve_diffusion(A[:,3], bval_hom, kdiff_A, Dxx, Dyy, GC_nmann)

    # exponentiate log terms
    n[:] = exp.(lnn[:,3])
    Te[:] = exp.(lnTe[:,3])
    Ti[:] = exp.(lnTi[:,3])
    Pe[:] = n .* Te
    Pi[:] = n .* Ti

    # solve vorticity and helmholtz equations
    lnn_x = Dx * lnn[:,3]
    lnn_y = Dy * lnn[:,3]
    Pi_xx = Dxx * Pi
    Pi_yy = Dyy * Pi
    phi = solve_vorticity_eqn(phi, phi_b, w[:,3], n, lnn_x, lnn_y, Pi_xx, Pi_yy, ad, Dx, Dy, Dxx, Dyy, GC_dchlt)
    psi = solve_helmholtz_eqn(psi, A[:,3], de2, Dxx, Dyy, GC_dchlt)

    # compute current
    j[:] = (Dxx + Dyy) * psi
    jn[:] = j ./ n

    # shift indices
    lnn[:,1:2] = lnn[:,2:3]
    lnTe[:,1:2] = lnTe[:,2:3]
    lnTi[:,1:2] = lnTi[:,2:3]
    u[:,1:2] = u[:,2:3]
    w[:,1:2] = w[:,2:3]
    A[:,1:2] = A[:,2:3]
end


# integrate N timesteps without saving
function leapfrog!(N, args...)

    for t=1:N
        leapfrog!(args...)
    end
end
