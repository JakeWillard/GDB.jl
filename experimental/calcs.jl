
struct Workspace

    GRID :: Grid
    Dx :: SparseMatrixCSC
    Dy :: SparseMatrixCSC
    Dxy :: SparseMatrixCSC
    Dxx :: SparseMatrixCSC
    Dyy :: SparseMatrixCSC
    Dxxx :: SparseMatrixCSC
    Dyyy :: SparseMatrixCSC
    Dxxy :: SparseMatrixCSC
    Dxyy :: SparseMatrixCSC
    Ds :: SparseMatrixCSC
    Dss :: SparseMatrixCSC
    DIFF_lnn :: LinearLeftHandSide
    DIFF_lnTe :: LinearLeftHandSide
    DIFF_lnTi :: LinearLeftHandSide
    DIFF_u :: LinearLeftHandSide
    DIFF_w :: LinearLeftHandSide
    DIFF_A :: LinearLeftHandSide
    HHOLTZ :: LinearLeftHandSide
    P0 :: SparseMatrixCSC
    P1 :: SparseMatrixCSC
    P2 :: SparseMatrixCSC
    P3 :: SparseMatrixCSC
    R1 :: SparseMatrixCSC
    R2 :: SparseMatrixCSC
    R3 :: SparseMatrixCSC
    LAM :: SparseMatrixCSC
    DCHLT1 :: SparseMatrixCSC
    DCHLT2 :: SparseMatrixCSC
    DCHLT3 :: SparseMatrixCSC
    NMANN1 :: SparseMatrixCSC
    NMANN2 :: SparseMatrixCSC
    NMANN3 :: SparseMatrixCSC
    FLXAVG :: SparseMatrixCSC
    TRGT :: Vector{Float64}
    Sn :: Vector{Float64}
    STe :: Vector{Float64}
    STi :: Vector{Float64}
    params :: Vector{Float64}
    dt :: Float64
    N_subcycle :: Int64

end


function lnn_forward_euler(lnn, lnn_t, wrk::Workspace)

    return wrk.LAM*(lnn + wrk.dt*(wrk.P0*lnn_t + wrk.P1*wrk.Sn))
end


function lnTe_forward_euler(lnTe, lnTe_t, wrk::Workspace)

    return wrk.LAM*(lnTe + wrk.dt*(wrk.P0*lnTe_t + wrk.P1*wrk.STe))
end


function lnTi_forward_euler(lnTi, lnTi_t, wrk::Workspace)

    return wrk.LAM*(lnTi + wrk.dt*(wrk.P0*lnTi_t + wrk.P1*wrk.STi))
end


function u_forward_euler(u, u_t, Te, Ti, wrk::Workspace)

    cs = wrk.TRGT .* sqrt.(Te + Ti)
    return wrk.LAM*(u + wrk.dt*(wrk.P0*u_t + wrk.P3*cs))
end


function w_forward_euler(w, w_t, wrk::Workspace)

    return wrk.LAM*(w + wrk.dt*wrk.P0*w_t)
end


function A_forward_euler(A, A_t, Te, Ti, n, phi, ev, ad, de2, wrk::Workspace)

    Abohm = wrk.TRGT .* Abohm_def(Te, Ti, n, phi, ev, ad, de2)

    return wrk.LAM*(A + wrk.dt*(wrk.P0*A_T + wrk.P3*Abohm))
end


function diffusion_lnn(lnn, wrk::Workspace)

    rhs = wrk.P0 * lnn
    return jacobi_preconditioned_gmres(wrk.DIFF_lnn, lnn, rhs)
end


function diffusion_lnTe(lnTe, wrk::Workspace)

    rhs = wrk.P0 * lnTe
    return jacobi_preconditioned_gmres(wrk.DIFF_lnTe, lnTe, rhs)
end


function diffusion_lnTi(lnTi, wrk::Workspace)

    rhs = wrk.P0 * lnTi
    return jacobi_preconditioned_gmres(wrk.DIFF_lnTi, lnTi, rhs)
end


function diffusion_u(u, Te, Ti, wrk::Workspace)

    cs = wrk.TRGT .* sqrt.(Te + Ti)
    rhs = wrk.P0 * u + wrk.DCHLT3*cs
    return jacobi_preconditioned_gmres(wrk.DIFF_u, u, rhs)
end


function diffusion_w(w, wrk::Workspace)

    rhs = wrk.P0 * w
    return jacobi_preconditioned_gmres(wrk.DIFF_w, w, rhs)
end


function diffusion_A(A, Te, Ti, n, phi, ev, ad, de2, wrk::Workspace)

    Abohm = wrk.TRGT .* Abohm_def(Te, Ti, n, phi, ev, ad, de2)
    rhs = wrk.P0 * A + wrk.DCHLT3*Abohm
    return jacobi_preconditioned_gmres(wrk.DIFF_A, A, rhs)
end


function vorticity_eqn(phi, w, Pi_xx, Pi_yy, Te, n, lnn_x, lnn_y, phi_b, ad, wrk::Workspace)

    M = Diagonal(n) * (wrk.Dxx + wrk.Dyy + Diagonal(lnn_x)*wrk.Dx + Diagonal(lnn_y)*wrk.Dy)
    LHS = LinearLeftHandSide(wrk.P0*M + wrk.DCHLT1 + wrk.DCHLT2 + wrk.DCHLT3)

    lambda = 2.695
    rhs0 = w - ad*(Pi_xx + Pi_yy)
    rhs = wrk.P0*rhs0 + wrk.DCHLT1*phi_b + lambda*(wrk.DCHLT2*Te + wrk.DCHLT3*Te)

    return jacobi_preconditioned_gmres(LHS, phi, rhs)
end


function helmholtz_eqn(psi, wrk::Workspace)

    rhs = wrk.P0 * psi
    return jacobi_preconditioned(wrk.HHOLTZ, psi, rhs)
end


function leapfrog!(lnn, lnTe, lnTi, u, w, A, phi, psi, n, Te, Ti, Pe, Pi, j, jn, wrk::Workspace)

    # unpack physical parameters
    am, ad, ki, ke, er, eg, ev, de2, eta = wrk.params

    # subcycle electron thermal conduction terms
    lnTe0 = lnTe[:,2]
    for _=1:wrk.N_subcycle

        # derivatives
        lnTe_x = wrk.Dx*lnTe[:,2]
        lnTe_y = wrk.Dy*lnTe[:,2]
        lnTe_xy = wrk.Dxy*lnTe[:,2]
        lnTe_xx = wrk.Dxx*lnTe[:,2]
        lnTe_yy = wrk.Dyy*lnTe[:,2]
        lnTe_s = partial_s(lnTe[:,2], lnTe_x, lnTe_y, psi_x, psi_y, wrk.Ds, am)
        lnTe_ss = partial_ss(lnTe[:,2], lnTe_x, lnTe_y, lnTe_xy, lnTe_xx, lnTe_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, wrk.Dss, wrk.Ds, wrk.Dx, wrk.Dy, am)

        lnTe_t = lnT_partial_thermal(lnTe_s, lnTe_ss, Te, n, ke)
        lnTe[:,3] = wrk.LAM *(0.5*(lnTe[:,1] + lnTe[:,2]) + wrk.dt*(wrk.P0*lnTe_t + wrk.P1*wrk.STe)/wrk.N_subcycle)

        # derivatives again
        lnTe_x = wrk.Dx*lnTe[:,2]
        lnTe_y = wrk.Dy*lnTe[:,2]
        lnTe_xy = wrk.Dxy*lnTe[:,2]
        lnTe_xx = wrk.Dxx*lnTe[:,2]
        lnTe_yy = wrk.Dyy*lnTe[:,2]
        lnTe_s = partial_s(lnTe[:,2], lnTe_x, lnTe_y, psi_x, psi_y, wrk.Ds, am)
        lnTe_ss = partial_ss(lnTe[:,2], lnTe_x, lnTe_y, lnTe_xy, lnTe_xx, lnTe_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, wrk.Dss, wrk.Ds, wrk.Dx, wrk.Dy, am)

        lnTe_t = lnT_partial_thermal(lnTe_s, lnTe_ss, Te, n, ke)
        lnTe[:,3] = wrk.LAM *(lnTe[:,2] + wrk.dt*(wrk.P0*lnTe_t + wrk.P1*wrk.STe)/wrk.N_subcycle)

        # shift indices
        lnTe[:,1:2] = lnTe[:,2:3]
    end
    lnTe[:,2] = 0.5*(lnTe0 + lnTe[:,2])

    # subcycle ion thermal conduction terms
    lnTi0 = lnTi[:,2]
    for _=1:wrk.N_subcycle

        # derivatives
        lnTi_x = wrk.Dx*lnTi[:,2]
        lnTi_y = wrk.Dy*lnTi[:,2]
        lnTi_xy = wrk.Dxy*lnTi[:,2]
        lnTi_xx = wrk.Dxx*lnTi[:,2]
        lnTi_yy = wrk.Dyy*lnTi[:,2]
        lnTi_s = partial_s(lnTi[:,2], lnTi_x, lnTi_y, psi_x, psi_y, wrk.Ds, am)
        lnTi_ss = partial_ss(lnTi[:,2], lnTi_x, lnTi_y, lnTi_xy, lnTi_xx, lnTi_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, wrk.Dss, wrk.Ds, wrk.Dx, wrk.Dy, am)

        lnTi_t = lnT_partial_thermal(lnTi_s, lnTi_ss, Te, n, ke)
        lnTi[:,3] = wrk.LAM *(0.5*(lnTi[:,1] + lnTi[:,2]) + wrk.dt*(wrk.P0*lnTi_t + wrk.P1*wrk.STe)/wrk.N_subcycle)

        # derivatives again
        lnTi_x = wrk.Dx*lnTi[:,2]
        lnTi_y = wrk.Dy*lnTi[:,2]
        lnTi_xy = wrk.Dxy*lnTi[:,2]
        lnTi_xx = wrk.Dxx*lnTi[:,2]
        lnTi_yy = wrk.Dyy*lnTi[:,2]
        lnTi_s = partial_s(lnTi[:,2], lnTi_x, lnTi_y, psi_x, psi_y, wrk.Ds, am)
        lnTi_ss = partial_ss(lnTi[:,2], lnTi_x, lnTi_y, lnTi_xy, lnTi_xx, lnTi_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, wrk.Dss, wrk.Ds, wrk.Dx, wrk.Dy, am)

        lnTi_t = lnT_partial_thermal(lnTi_s, lnTi_ss, Te, n, ke)
        lnTi[:,3] = wrk.LAM *(lnTi[:,2] + wrk.dt*(wrk.P0*lnTi_t + wrk.P1*wrk.STe)/wrk.N_subcycle)

        # shift indices
        lnTi[:,1:2] = lnTi[:,2:3]
    end
    lnTi[:,2] = 0.5*(lnTi0 + lnTi[:,2])

    # take derivatives
    psi_x = wrk.Dx*psi
    psi_y = wrk.Dy*psi
    psi_xy = wrk.Dxy*psi
    psi_xx = wrk.Dxx*psi
    psi_yy = wrk.Dyy*psi

    phi_x = wrk.Dx*phi
    phi_y = wrk.Dy*phi
    phi_xy = wrk.Dxy*phi
    phi_xx = wrk.Dxx*phi
    phi_yy = wrk.Dyy*phi
    phi_xxx = wrk.Dxxx*phi
    phi_yyy = wrk.Dyyy*phi
    phi_xxy = wrk.Dxxy*phi
    phi_xyy = wrk.Dxyy*phi
    phi_c = -2*phi_y[:]
    phi_s = partial_s(phi, phi_x, phi_y, psi_x, psi_y, wrk.Ds, am)

    Pe_x = wrk.Dx*Pe
    Pe_y = wrk.Dy*Pe
    Pe_xy = wrk.Dxy*Pe
    Pe_xx = wrk.Dxx*Pe
    Pe_yy = wrk.Dyy*Pe
    Pe_xxx = wrk.Dxxx*Pe
    Pe_yyy = wrk.Dyyy*Pe
    Pe_xxy = wrk.Dxxy*Pe
    Pe_xyy = wrk.Dxyy*Pe
    Pe_c = -2*Pe_y[:]
    Pe_s = partial_s(Pe, Pe_x, Pe_y, psi_x, psi_y, wrk.Ds, am)

    Pi_x = wrk.Dx*Pi
    Pi_y = wrk.Dy*Pi
    Pi_c = -2*Pi_y[:]
    Pi_s = partial_s(Pi, Pi_x, Pi_y, psi_x, psi_y, wrk.Ds, am)

    lnn_x = wrk.Dx*lnn[:,2]
    lnn_y = wrk.Dy*lnn[:,2]
    lnn_xy = wrk.Dxy*lnn[:,2]
    lnn_xx = wrk.Dxx*lnn[:,2]
    lnn_yy = wrk.Dyy*lnn[:,2]
    lnn_s = partial_s(lnn[:,2], lnn_x, lnn_y, psi_x, psi_y, wrk.Ds, am)

    lnTe_x = wrk.Dx*lnTe[:,2]
    lnTe_y = wrk.Dy*lnTe[:,2]
    lnTe_c = -2*lnTe_y[:]
    lnTe_s = partial_s(lnTe[:,2], lnTe_x, lnTe_y, psi_x, psi_y, wrk.Ds, am)

    lnTi_x = wrk.Dx*lnTi[:,2]
    lnTi_y = wrk.Dy*lnTi[:,2]
    lnTi_c = -2*lnTi_y[:]
    lnTi_s = partial_s(lnTi[:,2], lnTi_x, lnTi_y, psi_x, psi_y, wrk.Ds, am)

    u_x = wrk.Dx*u[:,2]
    u_y = wrk.Dy*u[:,2]
    u_c = -2*u_y[:]
    u_s = partial_s(u[:,2], u_x, u_y, psi_x, psi_y, wrk.Ds, am)

    j_x = wrk.Dx*j
    j_y = wrk.Dy*j
    j_s = partial_s(j, j_x, j_y, psi_x, psi_y, wrk.Ds, am)

    jn_x = wrk.Dx*jn
    jn_y = wrk.Dy*jn
    jn_s = partial_s(jn, jn_x, jn_y, psi_x, psi_y, wrk.Ds, am)

    G = G_def(Ti, phi_c, Pi_c, u_s, er, ev, ad)
    G_x = wrk.Dx*G
    G_y = wrk.Dy*G
    G_c = -2*G_y[:]

    # evaluate total derivatives
    lnn_Dt = lnn_total_ion_derivative(phi_c, Pe_c, n[:,2], j_s, u_s, er, ev, ad)
    lnTe_Dt = lnTe_total_electron_derivative_adiabatic(lnn_Dt, j, lnn_s, n[:,2], Te, lnTe_c, er, ad)
    lnTi_Dt = lnTi_total_ion_derivative_adiabatic(lnn_Dt, Ti, lnTi_c, er, ad)
    u_Dt = u_total_ion_derivative(Pe_s, Pi_s, G_s, n[:,2], Ti, u_c, ev, eg, er, ad)

    # evaluate partial derivatives
    lnn_t = general_partial_derivative(lnn_Dt, lnn_x, lnn_y, lnn_s, phi_x, phi_y, u[:,2])
    lnTe_t = general_partial_derivative(lnTe_Dt, lnTe_x, lnTe_y, lnTe_s, phi_x, phi_y, u[:,2] - jn / ev)
    lnTi_t = general_partial_derivative(lnTi_Dt, lnTi_x, lnTi_y, lnTi_s, phi_x, phi_y, u[:,2])
    u_t = general_partial_derivative(u_Dt, u_x, u_y, u_s, phi_x, phi_y, u[:,2])
    w_t = w_partial_t(n[:,2], lnn_x, lnn_y, lnn_xy, lnn_xx, lnn_yy,
                         Pe_x, Pe_y, Pe_xy, Pe_xx, Pe_yy, Pe_xxx, Pe_yyy, Pe_xxy, Pe_xyy, Pe_c,
                         phi_x, phi_y, phi_xy, phi_xx, phi_yy, phi_xxx, phi_yyy, phi_xxy, phi_xyy,
                         G_c, Pi_c, j_s, ad, eg)
    A_t = A_partial_t(phi_x, phi_y, phi_s, j, jn, jn_x, jn_y, jn_s, u[:,2], n[:,2], Te, de2, ev, ad, am, eta)

    # compute predictor step
    lnn[:,3] = lnn_forward_euler(0.5*(lnn[:,1] + lnn[:,2]), lnn_t, wrk)
    lnTe[:,3] = lnTe_forward_euler(0.5*(lnTe[:,1] + lnTe[:,2]), lnTe_t, wrk)
    lnTi[:,3] = lnTi_forward_euler(0.5*(lnTi[:,1] + lnTi[:,2]), lnTi_t, wrk)
    u[:,3] = u_forward_euler(0.5*(u[:,1] + u[:,2]), u_t, Te, Ti, wrk)
    w[:,3] = w_forward_euler(0.5*(w[:,1] + w[:,2]), w_t, wrk)
    A[:,3] = A_forward_euler(0.5*(A[:,1] + A[:,2]), A_t, Te, Ti, n[:,2], phi, ev, ad, de2, wrk)

    # compute inner boundary value for phi
    phi_b = dot(wrk.FLXAVG, phi)

    # solve linear problems
    lnn[:,3] = diffusion_lnn(lnn[:,3], wrk)
    lnTe[:,3] = diffusion_lnTe(lnTe[:,3], wrk)
    lnTi[:,3] = diffusion_lnTi(lnTi[:,3], wrk)
    u[:,3] = diffusion_u(u[:,3])
    w[:,3] = diffusion_w(w[:,3], wrk)
    phi = vorticity_eqn(phi, w[:,3], Pi_xx, Pi_yy, Te, n[:,3], lnn_x, lnn_y, phi_b, ad, wrk)
    psi = helmholtz_eqn(psi, wrk)

    n = exp.(lnn[:,3])
    Te = exp.(lnTe[:,3])
    Ti = exp.(lnTi[:,3])
    Pe = n .* Te
    Pi = n .* Ti
    j = (wrk.Dxx + wrk.Dyy) * psi
    jn = j ./ n

    # take derivatives
    psi_x = wrk.Dx*psi
    psi_y = wrk.Dy*psi
    psi_xy = wrk.Dxy*psi
    psi_xx = wrk.Dxx*psi
    psi_yy = wrk.Dyy*psi

    phi_x = wrk.Dx*phi
    phi_y = wrk.Dy*phi
    phi_xy = wrk.Dxy*phi
    phi_xx = wrk.Dxx*phi
    phi_yy = wrk.Dyy*phi
    phi_xxx = wrk.Dxxx*phi
    phi_yyy = wrk.Dyyy*phi
    phi_xxy = wrk.Dxxy*phi
    phi_xyy = wrk.Dxyy*phi
    phi_c = -2*phi_y[:]
    phi_s = partial_s(phi, phi_x, phi_y, psi_x, psi_y, wrk.Ds, am)

    Pe_x = wrk.Dx*Pe
    Pe_y = wrk.Dy*Pe
    Pe_xy = wrk.Dxy*Pe
    Pe_xx = wrk.Dxx*Pe
    Pe_yy = wrk.Dyy*Pe
    Pe_xxx = wrk.Dxxx*Pe
    Pe_yyy = wrk.Dyyy*Pe
    Pe_xxy = wrk.Dxxy*Pe
    Pe_xyy = wrk.Dxyy*Pe
    Pe_c = -2*Pe_y[:]
    Pe_s = partial_s(Pe, Pe_x, Pe_y, psi_x, psi_y, wrk.Ds, am)

    Pi_x = wrk.Dx*Pi
    Pi_y = wrk.Dy*Pi
    Pi_c = -2*Pi_y[:]
    Pi_s = partial_s(Pi, Pi_x, Pi_y, psi_x, psi_y, wrk.Ds, am)

    lnn_x = wrk.Dx*lnn[:,3]
    lnn_y = wrk.Dy*lnn[:,3]
    lnn_xy = wrk.Dxy*lnn[:,3]
    lnn_xx = wrk.Dxx*lnn[:,3]
    lnn_yy = wrk.Dyy*lnn[:,3]
    lnn_s = partial_s(lnn[:,3], lnn_x, lnn_y, psi_x, psi_y, wrk.Ds, am)

    lnTe_x = wrk.Dx*lnTe[:,3]
    lnTe_y = wrk.Dy*lnTe[:,3]
    lnTe_c = -2*lnTe_y[:]
    lnTe_s = partial_s(lnTe[:,3], lnTe_x, lnTe_y, psi_x, psi_y, wrk.Ds, am)

    lnTi_x = wrk.Dx*lnTi[:,3]
    lnTi_y = wrk.Dy*lnTi[:,3]
    lnTi_c = -2*lnTi_y[:]
    lnTi_s = partial_s(lnTi[:,3], lnTi_x, lnTi_y, psi_x, psi_y, wrk.Ds, am)

    u_x = wrk.Dx*u[:,3]
    u_y = wrk.Dy*u[:,3]
    u_c = -2*u_y[:]
    u_s = partial_s(u[:,3], u_x, u_y, psi_x, psi_y, wrk.Ds, am)

    j_x = wrk.Dx*j
    j_y = wrk.Dy*j
    j_s = partial_s(j, j_x, j_y, psi_x, psi_y, wrk.Ds, am)

    jn_x = wrk.Dx*jn
    jn_y = wrk.Dy*jn
    jn_s = partial_s(jn, jn_x, jn_y, psi_x, psi_y, wrk.Ds, am)

    G = G_def(Ti, phi_c, Pi_c, u_s, er, ev, ad)
    G_x = wrk.Dx*G
    G_y = wrk.Dy*G
    G_c = -2*G_y[:]

    # evaluate total derivatives
    lnn_Dt = lnn_total_ion_derivative(phi_c, Pe_c, n[:,3], j_s, u_s, er, ev, ad)
    lnTe_Dt = lnTe_total_electron_derivative_adiabatic(lnn_Dt, j, lnn_s, n[:,3], Te, lnTe_c, er, ad)
    lnTi_Dt = lnTi_total_ion_derivative_adiabatic(lnn_Dt, Ti, lnTi_c, er, ad)
    u_Dt = u_total_ion_derivative(Pe_s, Pi_s, G_s, n[:,3], Ti, u_c, ev, eg, er, ad)

    # evaluate partial derivatives
    lnn_t = general_partial_derivative(lnn_Dt, lnn_x, lnn_y, lnn_s, phi_x, phi_y, u[:,3])
    lnTe_t = general_partial_derivative(lnTe_Dt, lnTe_x, lnTe_y, lnTe_s, phi_x, phi_y, u[:,3] - jn / ev)
    lnTi_t = general_partial_derivative(lnTi_Dt, lnTi_x, lnTi_y, lnTi_s, phi_x, phi_y, u[:,3])
    u_t = general_partial_derivative(u_Dt, u_x, u_y, u_s, phi_x, phi_y, u[:,3])
    w_t = w_partial_t(n[:,3], lnn_x, lnn_y, lnn_xy, lnn_xx, lnn_yy,
                         Pe_x, Pe_y, Pe_xy, Pe_xx, Pe_yy, Pe_xxx, Pe_yyy, Pe_xxy, Pe_xyy, Pe_c,
                         phi_x, phi_y, phi_xy, phi_xx, phi_yy, phi_xxx, phi_yyy, phi_xxy, phi_xyy,
                         G_c, Pi_c, j_s, ad, eg)
    A_t = A_partial_t(phi_x, phi_y, phi_s, j, jn, jn_x, jn_y, jn_s, u[:,3], n[:,3], Te, de2, ev, ad, am, eta)

    # compute corrector step
    lnn[:,3] = lnn_forward_euler(lnn[:,2], lnn_t, wrk)
    lnTe[:,3] = lnTe_forward_euler(lnTe[:,2], lnTe_t, wrk)
    lnTi[:,3] = lnTi_forward_euler(lnTi[:,2], lnTi_t, wrk)
    u[:,3] = u_forward_euler(u[:,2], u_t, Te, Ti, wrk)
    w[:,3] = w_forward_euler(w[:,2], w_t, wrk)
    A[:,3] = A_forward_euler(A[:,2], A_t, Te, Ti, n, phi, ev, ad, de2, wrk)

    # compute inner boundary value for phi
    phi_b = dot(wrk.FLXAVG, phi)

    # solve linear problems
    lnn[:,3] = diffusion_lnn(lnn[:,3], wrk)
    lnTe[:,3] = diffusion_lnTe(lnTe[:,3], wrk)
    lnTi[:,3] = diffusion_lnTi(lnTi[:,3], wrk)
    u[:,3] = diffusion_u(u[:,3])
    w[:,3] = diffusion_w(w[:,3], wrk)
    phi = vorticity_eqn(phi, w[:,3], Pi_xx, Pi_yy, Te, n[:,3], lnn_x, lnn_y, phi_b, ad, wrk)
    psi = helmholtz_eqn(psi, wrk)

    # shift indices
    lnn[:,1:2] = lnn[:,2:3]
    lnTe[:,1:2] = lnTe[:,2:3]
    lnTi[:,1:2] = lnTi[:,2:3]
    u[:,1:2] = u[:,2:3]
    w[:,1:2] = w[:,2:3]
    A[:,1:2] = A[:,2:3]

    n = exp.(lnn[:,3])
    Te = exp.(lnTe[:,3])
    Ti = exp.(lnTi[:,3])
    Pe = n .* Te
    Pi = n .* Ti
    j = (wrk.Dxx + wrk.Dyy) * psi
    jn = j ./ n
end
