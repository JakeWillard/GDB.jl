

function lnn_forward_euler(lnn, lnn_t, a::Assets)

    return a.LAM*(lnn + a.dt*(a.P0*lnn_t + a.P1*a.Sn))
end


function lnTe_forward_euler(lnTe, lnTe_t, a::Assets)

    return a.LAM*(lnTe + a.dt*(a.P0*lnTe_t + a.P1*a.STe))
end


function lnTi_forward_euler(lnTi, lnTi_t, a::Assets)

    return a.LAM*(lnTi + a.dt*(a.P0*lnTi_t + a.P1*a.STi))
end


function u_forward_euler(u, u_t, Te, Ti, a::Assets)

    cs = a.TRGT .* sqrt.(Te + Ti)
    return a.LAM*(u + a.dt*(a.P0*u_t + a.P3*cs))
end


function w_forward_euler(w, w_t, a::Assets)

    return a.LAM*(w + a.dt*a.P0*w_t)
end


function A_forward_euler(A, A_t, Te, Ti, n, phi, er, ev, ad, de2, a::Assets)

    Abohm = a.TRGT .* Abohm_def(Te, Ti, n, phi, er, ev, ad, de2)

    return a.LAM*(A + a.dt*(a.P0*A_t))# + a.P3*Abohm))
end


function diffusion_lnn(lnn, a::Assets)

    rhs = a.P0 * lnn
    return a.DIFF_lnn.A \ rhs
    # return parallel_jacobi_preconditioned_gmres(a.DIFF_lnn, lnn, rhs, 1, a.GRID.Nz)
end


function diffusion_lnTe(lnTe, a::Assets)

    rhs = a.P0 * lnTe
    return a.DIFF_lnTe.A \ rhs
    # return parallel_jacobi_preconditioned_gmres(a.DIFF_lnTe, lnTe, rhs, 1, a.GRID.Nz)
end


function diffusion_lnTi(lnTi, a::Assets)

    rhs = a.P0 * lnTi
    return a.DIFF_lnTi.A \ rhs
    # return parallel_jacobi_preconditioned_gmres(a.DIFF_lnTi, lnTi, rhs, 1, a.GRID.Nz)
end


function diffusion_u(u, Te, Ti, a::Assets)

    cs = a.TRGT .* sqrt.(Te + Ti)
    rhs = a.P0 * u + a.DCHLT3*cs
    return a.DIFF_u.A \ rhs
    # return parallel_jacobi_preconditioned_gmres(a.DIFF_u, u, rhs, 1, a.GRID.Nz)
end


function diffusion_w(w, a::Assets)

    rhs = a.P0 * w
    return a.DIFF_w.A \ rhs
    # return parallel_jacobi_preconditioned_gmres(a.DIFF_w, w, rhs, 1, a.GRID.Nz)
end


function diffusion_A(A, Te, Ti, n, phi, er, ev, ad, de2, a::Assets)

    Abohm = a.TRGT .* Abohm_def(Te, Ti, n, phi, er, ev, ad, de2)
    rhs = a.P0 * A # + a.DCHLT3*Abohm
    return a.DIFF_A.A \ rhs
    # return parallel_jacobi_preconditioned_gmres(a.DIFF_A, A, rhs, 1, a.GRID.Nz)
end


function vorticity_eqn(phi, w, Pi_xx, Pi_yy, Te, n, lnn_x, lnn_y, phi_b, ad, a::Assets)

    M = (a.Dxx + a.Dyy)# + Diagonal(lnn_x)*a.Dx + Diagonal(lnn_y)*a.Dy)
    LHS = LinearLeftHandSide(a.P0*M + a.DCHLT1 + a.DCHLT2 + a.DCHLT3, 2/3)

    lambda = 2.695
    rhs0 = (w - ad*(Pi_xx + Pi_yy)) ./ n
    rhs = a.P0*rhs0 + a.DCHLT1*phi_b + lambda*(a.DCHLT2*Te + a.DCHLT3*Te)
    return LHS.A \ rhs
    # return parallel_jacobi_preconditioned_gmres(LHS, phi, rhs, 1, a.GRID.Nz; m=100)
end


function helmholtz_eqn(A, a::Assets)

    rhs = a.P0 * A
    return a.HHOLTZ.A \ rhs
    # return parallel_jacobi_preconditioned_gmres(a.HHOLTZ, psi, rhs, 1, a.GRID.Nz)
end


function leapfrog!(lnn, lnTe, lnTi, u, w, A, phi, psi, n, Te, Ti, Pe, Pi, j, jn, a::Assets)

    # unpack physical parameters
    am, ad, ki, ke, er, eg, ev, de2, eta = a.params

    # take derivatives of psi
    psi_x = a.Dx*psi
    psi_y = a.Dy*psi
    psi_xy = a.Dxy*psi
    psi_xx = a.Dxx*psi
    psi_yy = a.Dyy*psi

    # subcycle electron thermal conduction terms
    lnTe0 = lnTe[:,2]
    for _=1:a.N_subcycle

        # derivatives
        lnTe_x = a.Dx*lnTe[:,2]
        lnTe_y = a.Dy*lnTe[:,2]
        lnTe_xy = a.Dxy*lnTe[:,2]
        lnTe_xx = a.Dxx*lnTe[:,2]
        lnTe_yy = a.Dyy*lnTe[:,2]
        lnTe_s = partial_s(lnTe[:,2], lnTe_x, lnTe_y, psi_x, psi_y, a.Ds, am)
        lnTe_ss = partial_ss(lnTe[:,2], lnTe_s, lnTe_x, lnTe_y, lnTe_xy, lnTe_xx, lnTe_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, a.Dss, a.Ds, a.Dx, a.Dy, am)

        lnTe_t = lnT_partial_thermal(lnTe_s, lnTe_ss, Te, n, ke)
        lnTe[:,3] = a.LAM *(0.5*(lnTe[:,1] + lnTe[:,2]) + a.dt*(a.P0*lnTe_t + a.P1*a.STe)/a.N_subcycle)

        # derivatives again
        lnTe_x = a.Dx*lnTe[:,2]
        lnTe_y = a.Dy*lnTe[:,2]
        lnTe_xy = a.Dxy*lnTe[:,2]
        lnTe_xx = a.Dxx*lnTe[:,2]
        lnTe_yy = a.Dyy*lnTe[:,2]
        lnTe_s = partial_s(lnTe[:,2], lnTe_x, lnTe_y, psi_x, psi_y, a.Ds, am)
        lnTe_ss = partial_ss(lnTe[:,2], lnTe_s, lnTe_x, lnTe_y, lnTe_xy, lnTe_xx, lnTe_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, a.Dss, a.Ds, a.Dx, a.Dy, am)

        lnTe_t = lnT_partial_thermal(lnTe_s, lnTe_ss, Te, n, ke)
        lnTe[:,3] = a.LAM *(lnTe[:,2] + a.dt*(a.P0*lnTe_t + a.P1*a.STe)/a.N_subcycle)

        # shift indices
        lnTe[:,1:2] = lnTe[:,2:3]
    end
    lnTe[:,2] = 0.5*(lnTe0 + lnTe[:,2])

    # subcycle ion thermal conduction terms
    lnTi0 = lnTi[:,2]
    for _=1:a.N_subcycle

        # derivatives
        lnTi_x = a.Dx*lnTi[:,2]
        lnTi_y = a.Dy*lnTi[:,2]
        lnTi_xy = a.Dxy*lnTi[:,2]
        lnTi_xx = a.Dxx*lnTi[:,2]
        lnTi_yy = a.Dyy*lnTi[:,2]
        lnTi_s = partial_s(lnTi[:,2], lnTi_x, lnTi_y, psi_x, psi_y, a.Ds, am)
        lnTi_ss = partial_ss(lnTi[:,2], lnTi_s, lnTi_x, lnTi_y, lnTi_xy, lnTi_xx, lnTi_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, a.Dss, a.Ds, a.Dx, a.Dy, am)

        lnTi_t = lnT_partial_thermal(lnTi_s, lnTi_ss, Te, n, ke)
        lnTi[:,3] = a.LAM *(0.5*(lnTi[:,1] + lnTi[:,2]) + a.dt*(a.P0*lnTi_t + a.P1*a.STe)/a.N_subcycle)

        # derivatives again
        lnTi_x = a.Dx*lnTi[:,2]
        lnTi_y = a.Dy*lnTi[:,2]
        lnTi_xy = a.Dxy*lnTi[:,2]
        lnTi_xx = a.Dxx*lnTi[:,2]
        lnTi_yy = a.Dyy*lnTi[:,2]
        lnTi_s = partial_s(lnTi[:,2], lnTi_x, lnTi_y, psi_x, psi_y, a.Ds, am)
        lnTi_ss = partial_ss(lnTi[:,2], lnTi_s, lnTi_x, lnTi_y, lnTi_xy, lnTi_xx, lnTi_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, a.Dss, a.Ds, a.Dx, a.Dy, am)

        lnTi_t = lnT_partial_thermal(lnTi_s, lnTi_ss, Te, n, ke)
        lnTi[:,3] = a.LAM *(lnTi[:,2] + a.dt*(a.P0*lnTi_t + a.P1*a.STe)/a.N_subcycle)

        # shift indices
        lnTi[:,1:2] = lnTi[:,2:3]
    end
    lnTi[:,2] = 0.5*(lnTi0 + lnTi[:,2])

    # take other derivatives
    phi_x = a.Dx*phi
    phi_y = a.Dy*phi
    phi_xy = a.Dxy*phi
    phi_xx = a.Dxx*phi
    phi_yy = a.Dyy*phi
    phi_xxx = a.Dxxx*phi
    phi_yyy = a.Dyyy*phi
    phi_xxy = a.Dxxy*phi
    phi_xyy = a.Dxyy*phi
    phi_c = -2*phi_y[:]
    phi_s = partial_s(phi, phi_x, phi_y, psi_x, psi_y, a.Ds, am)

    Pe_x = a.Dx*Pe
    Pe_y = a.Dy*Pe
    Pe_xy = a.Dxy*Pe
    Pe_xx = a.Dxx*Pe
    Pe_yy = a.Dyy*Pe
    Pe_xxx = a.Dxxx*Pe
    Pe_yyy = a.Dyyy*Pe
    Pe_xxy = a.Dxxy*Pe
    Pe_xyy = a.Dxyy*Pe
    Pe_c = -2*Pe_y[:]
    Pe_s = partial_s(Pe, Pe_x, Pe_y, psi_x, psi_y, a.Ds, am)

    Pi_x = a.Dx*Pi
    Pi_y = a.Dy*Pi
    Pi_c = -2*Pi_y[:]
    Pi_s = partial_s(Pi, Pi_x, Pi_y, psi_x, psi_y, a.Ds, am)

    lnn_x = a.Dx*lnn[:,2]
    lnn_y = a.Dy*lnn[:,2]
    lnn_xy = a.Dxy*lnn[:,2]
    lnn_xx = a.Dxx*lnn[:,2]
    lnn_yy = a.Dyy*lnn[:,2]
    lnn_s = partial_s(lnn[:,2], lnn_x, lnn_y, psi_x, psi_y, a.Ds, am)

    lnTe_x = a.Dx*lnTe[:,2]
    lnTe_y = a.Dy*lnTe[:,2]
    lnTe_c = -2*lnTe_y[:]
    lnTe_s = partial_s(lnTe[:,2], lnTe_x, lnTe_y, psi_x, psi_y, a.Ds, am)

    lnTi_x = a.Dx*lnTi[:,2]
    lnTi_y = a.Dy*lnTi[:,2]
    lnTi_c = -2*lnTi_y[:]
    lnTi_s = partial_s(lnTi[:,2], lnTi_x, lnTi_y, psi_x, psi_y, a.Ds, am)

    u_x = a.Dx*u[:,2]
    u_y = a.Dy*u[:,2]
    u_c = -2*u_y[:]
    u_s = partial_s(u[:,2], u_x, u_y, psi_x, psi_y, a.Ds, am)

    j_x = a.Dx*j
    j_y = a.Dy*j
    j_s = partial_s(j, j_x, j_y, psi_x, psi_y, a.Ds, am)

    jn_x = a.Dx*jn
    jn_y = a.Dy*jn
    jn_s = partial_s(jn, jn_x, jn_y, psi_x, psi_y, a.Ds, am)

    G = G_def(Ti, phi_c, Pi_c, u_s, er, ev, ad)
    G_x = a.Dx*G
    G_y = a.Dy*G
    G_c = -2*G_y[:]
    G_s = partial_s(G, G_x, G_y, psi_x, psi_y, a.Ds, am)

    # evaluate total derivatives
    lnn_Dt = lnn_total_ion_derivative(phi_c, Pe_c, n, j_s, u_s, er, ev, ad)
    lnTe_Dt = lnTe_total_electron_derivative_adiabatic(lnn_Dt, j, lnn_s, n, Te, lnTe_c, er, ad)
    lnTi_Dt = lnTi_total_ion_derivative_adiabatic(lnn_Dt, Ti, lnTi_c, er, ad)
    u_Dt = u_total_ion_derivative(Pe_s, Pi_s, G_s, n, Ti, u_c, ev, eg, er, ad)

    # evaluate partial derivatives
    lnn_t = general_partial_derivative(lnn_Dt, lnn_x, lnn_y, lnn_s, phi_x, phi_y, u[:,2])
    lnTe_t = general_partial_derivative(lnTe_Dt, lnTe_x, lnTe_y, lnTe_s, phi_x, phi_y, u[:,2] - jn / ev)
    lnTi_t = general_partial_derivative(lnTi_Dt, lnTi_x, lnTi_y, lnTi_s, phi_x, phi_y, u[:,2])
    u_t = general_partial_derivative(u_Dt, u_x, u_y, u_s, phi_x, phi_y, u[:,2])
    w_t = w_partial_t(n, lnn_x, lnn_y, lnn_xy, lnn_xx, lnn_yy,
                         Pe_x, Pe_y, Pe_xy, Pe_xx, Pe_yy, Pe_xxx, Pe_yyy, Pe_xxy, Pe_xyy, Pe_c,
                         phi_x, phi_y, phi_xy, phi_xx, phi_yy, phi_xxx, phi_yyy, phi_xxy, phi_xyy,
                         G_c, Pi_c, j_s, ad, eg)
    A_t = A_partial_t(phi_x, phi_y, phi_s, j, jn, jn_x, jn_y, jn_s, Pe_s, u[:,2], n, Te, de2, ev, ad, am, eta)

    # compute predictor step
    lnn[:,3] = lnn_forward_euler(0.5*(lnn[:,1] + lnn[:,2]), lnn_t, a)
    lnTe[:,3] = lnTe_forward_euler(0.5*(lnTe[:,1] + lnTe[:,2]), lnTe_t, a)
    lnTi[:,3] = lnTi_forward_euler(0.5*(lnTi[:,1] + lnTi[:,2]), lnTi_t, a)
    u[:,3] = u_forward_euler(0.5*(u[:,1] + u[:,2]), u_t, Te, Ti, a)
    w[:,3] = w_forward_euler(0.5*(w[:,1] + w[:,2]), w_t, a)
    A[:,3] = A_forward_euler(0.5*(A[:,1] + A[:,2]), A_t, Te, Ti, n, phi, er, ev, ad, de2, a)

    # compute inner boundary value for phi
    phi_b = dot(a.FLXAVG, phi) * ones(a.GRID.Nk*a.GRID.Nz)

    # solve linear problems
    lnn[:,3] = diffusion_lnn(lnn[:,3], a)
    lnTe[:,3] = diffusion_lnTe(lnTe[:,3], a)
    lnTi[:,3] = diffusion_lnTi(lnTi[:,3], a)
    u[:,3] = diffusion_u(u[:,3], Te, Ti, a)
    w[:,3] = diffusion_w(w[:,3], a)
    phi[:] = vorticity_eqn(phi, w[:,3], Pi_xx, Pi_yy, Te, n, lnn_x, lnn_y, phi_b, ad, a)
    psi[:] = helmholtz_eqn(A[:,3], a)

    n[:] = exp.(lnn[:,3])
    Te[:] = exp.(lnTe[:,3])
    Ti[:] = exp.(lnTi[:,3])
    Pe[:] = n .* Te
    Pi[:] = n .* Ti
    j[:] = (a.Dxx + a.Dyy) * psi
    jn[:] = j ./ n

    # take derivatives
    psi_x = a.Dx*psi
    psi_y = a.Dy*psi
    psi_xy = a.Dxy*psi
    psi_xx = a.Dxx*psi
    psi_yy = a.Dyy*psi

    phi_x = a.Dx*phi
    phi_y = a.Dy*phi
    phi_xy = a.Dxy*phi
    phi_xx = a.Dxx*phi
    phi_yy = a.Dyy*phi
    phi_xxx = a.Dxxx*phi
    phi_yyy = a.Dyyy*phi
    phi_xxy = a.Dxxy*phi
    phi_xyy = a.Dxyy*phi
    phi_c = -2*phi_y[:]
    phi_s = partial_s(phi, phi_x, phi_y, psi_x, psi_y, a.Ds, am)

    Pe_x = a.Dx*Pe
    Pe_y = a.Dy*Pe
    Pe_xy = a.Dxy*Pe
    Pe_xx = a.Dxx*Pe
    Pe_yy = a.Dyy*Pe
    Pe_xxx = a.Dxxx*Pe
    Pe_yyy = a.Dyyy*Pe
    Pe_xxy = a.Dxxy*Pe
    Pe_xyy = a.Dxyy*Pe
    Pe_c = -2*Pe_y[:]
    Pe_s = partial_s(Pe, Pe_x, Pe_y, psi_x, psi_y, a.Ds, am)

    Pi_x = a.Dx*Pi
    Pi_y = a.Dy*Pi
    Pi_c = -2*Pi_y[:]
    Pi_s = partial_s(Pi, Pi_x, Pi_y, psi_x, psi_y, a.Ds, am)

    lnn_x = a.Dx*lnn[:,3]
    lnn_y = a.Dy*lnn[:,3]
    lnn_xy = a.Dxy*lnn[:,3]
    lnn_xx = a.Dxx*lnn[:,3]
    lnn_yy = a.Dyy*lnn[:,3]
    lnn_s = partial_s(lnn[:,3], lnn_x, lnn_y, psi_x, psi_y, a.Ds, am)

    lnTe_x = a.Dx*lnTe[:,3]
    lnTe_y = a.Dy*lnTe[:,3]
    lnTe_c = -2*lnTe_y[:]
    lnTe_s = partial_s(lnTe[:,3], lnTe_x, lnTe_y, psi_x, psi_y, a.Ds, am)

    lnTi_x = a.Dx*lnTi[:,3]
    lnTi_y = a.Dy*lnTi[:,3]
    lnTi_c = -2*lnTi_y[:]
    lnTi_s = partial_s(lnTi[:,3], lnTi_x, lnTi_y, psi_x, psi_y, a.Ds, am)

    u_x = a.Dx*u[:,3]
    u_y = a.Dy*u[:,3]
    u_c = -2*u_y[:]
    u_s = partial_s(u[:,3], u_x, u_y, psi_x, psi_y, a.Ds, am)

    j_x = a.Dx*j
    j_y = a.Dy*j
    j_s = partial_s(j, j_x, j_y, psi_x, psi_y, a.Ds, am)

    jn_x = a.Dx*jn
    jn_y = a.Dy*jn
    jn_s = partial_s(jn, jn_x, jn_y, psi_x, psi_y, a.Ds, am)

    G = G_def(Ti, phi_c, Pi_c, u_s, er, ev, ad)
    G_x = a.Dx*G
    G_y = a.Dy*G
    G_c = -2*G_y[:]
    G_s = partial_s(G, G_x, G_y, psi_x, psi_y, a.Ds, am)

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

    # compute corrector step
    lnn[:,3] = lnn_forward_euler(lnn[:,2], lnn_t, a)
    lnTe[:,3] = lnTe_forward_euler(lnTe[:,2], lnTe_t, a)
    lnTi[:,3] = lnTi_forward_euler(lnTi[:,2], lnTi_t, a)
    u[:,3] = u_forward_euler(u[:,2], u_t, Te, Ti, a)
    w[:,3] = w_forward_euler(w[:,2], w_t, a)
    A[:,3] = A_forward_euler(A[:,2], A_t, Te, Ti, n, phi, er, ev, ad, de2, a)

    # compute inner boundary value for phi
    phi_b = dot(a.FLXAVG, phi) * ones(a.GRID.Nk*a.GRID.Nz)

    # solve linear problems
    lnn[:,3] = diffusion_lnn(lnn[:,3], a)
    lnTe[:,3] = diffusion_lnTe(lnTe[:,3], a)
    lnTi[:,3] = diffusion_lnTi(lnTi[:,3], a)
    u[:,3] = diffusion_u(u[:,3], Te, Ti, a)
    w[:,3] = diffusion_w(w[:,3], a)
    phi[:] = vorticity_eqn(phi, w[:,3], Pi_xx, Pi_yy, Te, n, lnn_x, lnn_y, phi_b, ad, a)
    psi[:] = helmholtz_eqn(A[:,3], a)

    # shift indices
    lnn[:,1:2] = lnn[:,2:3]
    lnTe[:,1:2] = lnTe[:,2:3]
    lnTi[:,1:2] = lnTi[:,2:3]
    u[:,1:2] = u[:,2:3]
    w[:,1:2] = w[:,2:3]
    A[:,1:2] = A[:,2:3]

    n[:] = exp.(lnn[:,3])
    Te[:] = exp.(lnTe[:,3])
    Ti[:] = exp.(lnTi[:,3])
    Pe[:] = n .* Te
    Pi[:] = n .* Ti
    j[:] = (a.Dxx + a.Dyy) * psi
    jn[:] = j ./ n
end
