
struct TimestepConstants

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
    Sn :: Vector{Float64}
    STe :: Vector{Float64}
    STi :: Vector{Float64}
    dt :: Float64
    N_subcycle :: Int64

end


function lnn_forward_euler(lnn, lnn_t, tsc::TimestepConstants)

    return tsc.LAM*(lnn + tsc.dt*(tsc.P0*lnn_t + tsc.P1*tcs.Sn))
end


function lnTe_forward_euler(lnTe, lnTe_t, tsc::TimestepConstants)

    return tsc.LAM*(lnTe + tsc.dt*(tsc.P0*lnTe_t + tsc.P1*tcs.STe))
end


function lnTi_forward_euler(lnTi, lnTi_t, tsc::TimestepConstants)

    return tsc.LAM*(lnTi + tsc.dt*(tsc.P0*lnTi_t + tsc.P1*tcs.STi))
end


function u_forward_euler(u, u_t, Te, Ti, tsc::TimestepConstants)

    cs = sqrt.(Te + Ti)
    return tsc.LAM*(u + tcs.dt*(tsc.P0*u_t + tcs.P3*cs))
end


function w_forward_euler(w, w_t, tsc::TimestepConstants)

    return tcs.LAM*(w + tcs.dt*tcs.P0*w_t)
end


function A_forward_euler(A, A_t, Te, Ti, n, phi, ev, ad, de2, tsc::TimestepConstants)

    Abohm = Abohm_def(Te, Ti, n, phi, ev, ad, de2)

    return tcs.LAM*(A + tcs.dt*(tcs.P0*A_T + tcs.P3*Abohm))
end


function diffusion_lnn(lnn, tcs::TimestepConstants)

    rhs = tcs.P0 * lnn
    return jacobi_preconditiond_gmres(tcs.DIFF_lnn, lnn, rhs)
end


function diffusion_lnTe(lnTe, tcs::TimestepConstants)

    rhs = tcs.P0 * lnTe
    return jacobi_preconditiond_gmres(tcs.DIFF_lnTe, lnTe, rhs)
end


function diffusion_lnTi(lnTi, tcs::TimestepConstants)

    rhs = tcs.P0 * lnTi
    return jacobi_preconditiond_gmres(tcs.DIFF_lnTi, lnTi, rhs)
end


function diffusion_u(u, Te, Ti, tcs::TimestepConstants)

    cs = sqrt.(Te + Ti)
    rhs = tcs.P0 * u + tcs.DCHLT3*cs
    return jacobi_preconditiond_gmres(tcs.DIFF_u, u, rhs)
end


function diffusion_w(w, tcs::TimestepConstants)

    rhs = tcs.P0 * w
    return jacobi_preconditiond_gmres(tcs.DIFF_w, w, rhs)
end


function diffusion_A(A, Te, Ti, n, phi, ev, ad, de2, tcs::TimestepConstants)

    Abohm = Abohm_def(Te, Ti, n, phi, ev, ad, de2)
    rhs = tcs.P0 * A + tcs.DCHLT3*Abohm
    return jacobi_preconditiond_gmres(tcs.DIFF_A, A, rhs)
end


function vorticity_eqn(phi, w, Pi_xx, Pi_yy, Te, n, lnn_x, lnn_y, phi_b, ad, tcs::TimestepConstants)

    M = Diagonal(n) * (tcs.Dxx + tcs.Dyy + Diagonal(lnn_x)*tcs.Dx + Diagonal(lnn_y)*tcs.Dy)
    LHS = LinearLeftHandSide(tcs.P0*M + tcs.DCHLT1 + tcs.DCHLT2 + tcs.DCHLT3)

    lambda = 2.695
    rhs0 = w - ad*(Pi_xx + Pi_yy)
    rhs = tcs.P0*rhs0 + tcs.DCHLT1*phi_b + lambda*(tcs.DCHLT2*Te + tcs.DCHLT3*Te)

    return jacobi_preconditiond_gmres(LHS, phi, rhs)


function adiabatic_leapfrog!()

    psi_x = Dx*psi[:,2]
    psi_y = Dy*psi[:,2]
    psi_xy = Dxy*psi[:,2]
    psi_xx = Dxx*psi[:,2]
    psi_yy = Dyy*psi[:,2]

    phi_x = Dx*phi[:,2]
    phi_y = Dy*phi[:,2]
    phi_xy = Dxy*phi[:,2]
    phi_xx = Dxx*phi[:,2]
    phi_yy = Dyy*phi[:,2]
    phi_xxx = Dxxx*phi[:,2]
    phi_yyy = Dyyy*phi[:,2]
    phi_xxy = Dxxy*phi[:,2]
    phi_xyy = Dxyy*phi[:,2]
    phi_c = -2*phi_y[:]
    phi_s = partial_s(phi[:,2], phi_x, phi_y, psi_x, psi_y, Ds, am)

    Pe_x = Dx*Pe[:,2]
    Pe_y = Dy*Pe[:,2]
    Pe_xy = Dxy*Pe[:,2]
    Pe_xx = Dxx*Pe[:,2]
    Pe_yy = Dyy*Pe[:,2]
    Pe_xxx = Dxxx*Pe[:,2]
    Pe_yyy = Dyyy*Pe[:,2]
    Pe_xxy = Dxxy*Pe[:,2]
    Pe_xyy = Dxyy*Pe[:,2]
    Pe_c = -2*Pe_y[:]
    Pe_s = partial_s(Pe[:,2], Pe_x, Pe_y, psi_x, psi_y, Ds, am)

    Pi_x = Dx*Pi[:,2]
    Pi_y = Dy*Pi[:,2]
    Pi_c = -2*Pi_y[:]
    Pi_s = partial_s(Pi[:,2], Pi_x, Pi_y, psi_x, psi_y, Ds, am)

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

    j_x = Dx*j[:,2]
    j_y = Dy*j[:,2]
    j_s = partial_s(j[:,2], j_x, j_y, psi_x, psi_y, Ds, am)

    jn_x = Dx*jn
    jn_y = Dy*jn
    jn_s = partial_s(jn[:,2], jn_x, jn_y, psi_x, psi_y, Ds, am)

    G = G_def(Ti, phi_c, Pi_c, u_s, er, ev, ad)
    G_x = Dx*G
    G_y = Dy*G
    G_c = -2*G_y[:]

    # evaluate total derivatives
    lnn_Dt = lnn_total_ion_derivative(phi_c, Pe_c, n, j_s, u_s, er, ev, ad)
    lnTe_Dt = lnTe_total_electron_derivative_adiabatic(lnn_Dt, j, lnn_s, n, Te, lnTe_c, er, ad)
    lnTi_Dt = lnTi_total_ion_derivative_adiabatic(lnn_Dt, Ti, lnTi_c, er, ad)
    u_Dt = u_total_ion_derivative(Pe_s, Pi_s, n, G_s, n, Ti, u_c, ev, eg, er, ad)

    # evaluate partial derivatives
    lnn_t = general_partial_derivative(lnn_Dt, lnn_x, lnn_y, lnn_s, phi_x, phi_y, u)
    lnTe_t = general_partial_derivative(lnTe_Dt, lnTe_x, lnTe_y, lnTe_s, phi_x, phi_y, u - jn / ev)
    lnTi_t = general_partial_derivative(lnTi_Dt, lnTi_x, lnTi_y, lnTi_s, phi_x, phi_y, u)
    u_t = general_partial_derivative(u_Dt, u_x, u_y, u_s, phi_x, phi_y, u)
    w_t = w_partial_t(n, lnn_x, lnn_y, lnn_xy, lnn_xx, lnn_yy,
                         Pe_x, Pe_y, Pe_xy, Pe_xx, Pe_yy, Pe_xxx, Pe_yyy, Pe_xxy, Pe_xyy, Pe_c,
                         phi_x, phi_y, phi_xy, phi_xx, phi_yy, phi_xxx, phi_yyy, phi_xxy, phi_xyy,
                         G_c, Pi_c, j_s, ad, eg)
    A_t = A_partial_t(phi_x, phi_y, phi_s, j, jn, jn_x, jn_y, jn_s, u, n, Te, de2, ev, ad, am, eta)

    # compute predictor step
    lnn[:,3] = lnn_forward_euler(0.5*(lnn[:,1] + lnn[:,2]), lnn_t, Sn, P0, P1, LAM, dt)
    lnTe[:,3] = lnTe_forward_euler(0.5*(lnTe[:,1] + lnTe[:,2]), lnTe_t, STe, P0, P1, LAM, dt)
    lnTi[:,3] = lnTi_forward_euler(0.5*(lnTi[:,1] + lnTi[:,2]), lnTi_t, STi, P0, P1, LAM, dt)
    u[:,3] = u_forward_euler(0.5*(u[:,1] + u[:,2]), u_t, Te, Ti, P0, P3, LAM, dt)
    w[:,3] = w_forward_euler(0.5*(w[:,1] + w[:,2]), w_t, P0, LAM, dt)
    A[:,3] = A_forward_euler(0.5*(A[:,1] + A[:,2]), A_t, Te, Ti, n, phi, ev, ad, de2, P0, P3, LAM, dt)

    # solve linear problems
    lnn[:,3] = diffusion_lnn()
    lnTe[:,3] = diffusion_lnTe()
    lnTi[:,3] = diffusion_lnTi()
    u[:,3] = diffusion_u()
    w[:,3] = diffusion_w()
    phi = vorticity_eqn()
    psi = helmholtz_eqn()
    n = exp.(lnn[:,3])
    Te = exp.(lnTe[:,3])
    Ti = exp.(lnTi[:,3])
    j = (Dxx + Dyy) * psi
    jn = j ./ n 

    # take derivatives again
    # XXX

    # partial time derivatives
    lnn_t, lnTe_t, lnTi_t, u_t, w_t, A_t = all_partial_t()

    # compute corrector step
    lnn_forward_euler(lnn[:,2], lnn_t, Sn, P0, P1, LAM, dt)
    lnTe_forward_euler(lnTe[:,2], lnTe_t, STe, P0, P1, LAM, dt)
    lnTi_forward_euler(lnTi[:,2], lnTi_t, STi, P0, P1, LAM, dt)
    u_forward_euler(u[:,2], u_t, Te, Ti, P0, P3, LAM, dt)
    w_forward_euler(w[:,2], w_t, P0, LAM, dt)
    A_forward_euler(A[:,2], A_t, Te, Ti, n, phi, ev, ad, de2, P0, P3, LAM, dt)

    # solve linear problems
    # XXX
end
