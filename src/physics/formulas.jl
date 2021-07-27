

function partial_s(f, f_x, f_y, psi_x, psi_y, Ds, am)

    return Ds*f + am*(psi_x .* f_y - psi_y .* f_x)
end


function partial_ss()

    #XXX just apply partial_s twice for now to compute f_ss
end


function lnn_total_ion_derivative(phi_c, Pe_c, n, j_s, u_s, er, ev, ad)

    _a = er*phi_c - er*ad*Pe_c ./ n
    _b = er*ad*j_s ./ n - ev*u_s
    return _a + _b
end


function lnTe_total_electron_derivative_adiabatic(lnn_Dt, j, lnn_s, n, Te, lnTe_c, er, ad)

    _a = (2/3)*(lnn_Dt + er*ad*j.*lnn_s./n)
    _b = -(5/3)*er*ad*(Te .* lnTe_c)

    return _a + _b
end


function lnTi_total_ion_derivative_adiabatic(lnn_Dt, Ti, lnTi_c, er, ad)

    _a = (2/3)*lnn_Dt
    _b = (5/3)*er*ad*(Ti .* lnTi_c)

    return _a + _b
end


# NOTE: this one is written to be species-ambiguous since it applies to Te and Ti
function lnT_partial_thermal(lnT_s, lnT_ss, T, n, k)

    _a = (7/2)*lnT_s.^2 + lnT_ss
    return k*T.^(5/2) * _a ./ n
end


function u_total_ion_derivative(Pe_s, Pi_s, n, G_s, n, Ti, u_c, ev, eg, er, ad)

    _a = -ev*(Pe_s + Pi_s) ./ n
    _b = -4*ev*eg*G_s ./ n
    _c = er*ad*(Ti .* u_c)

    return +(_a, _b, _c)
end


function general_partial_derivative(f_Dt, f_x, f_y, f_s, phi_x, phi_y, v)

    return f_Dt - phi_x.*f_y + phi_y.*f_x - v*f_s
end


function w_partial_t(n, lnn_x, ) #XXX No dependence on lnn_xx, despite dependence on lnn_yy? Mistake?

    _a = n.*(lnn_x.*phi_x.*phi_xy + phi_xx.*phi_xy + phi_x.*phi_xxy)
    _b = -n.*(lnn_x.*phi_y*phi_xx + phi_xy.*phi_xx + phi_y.*phi_xxx)
    _c = n.*(lnn_y.*phi_x.*phi_yy + phi_xy.*phi_yy + phi_x.*phi_yyy)
    _d = -n.*(lnn_y.*phi_y.*phi_xy + phi_yy.*phi_xy + phi_y.*phi_xyy)
    _e = ad*(phi_xx.*Pi_xy + phi_x.*Pi_xxy - phi_xy.*Pi_xx - phi_y.*Pi_xxx)
    _f = ad*(phi_xx.*lnn_y.*Pi_x + phi_x.*lnn_xy.*Pi_x + phi_x.*lnn_y.*Pi_xx)
    _g = -ad*(phi_xy.*lnn_x.*Pi_x + phi_y.*lnn_xy.*Pi_x + phi_y.*lnn_x.*Pi_xx)
    _h = ad*(phi_xy.*Pi_yy + phi_x.*Pi_yyy - phi_yy.*Pi_xy - phi_y.*Pi_xyy)
    _i = ad*(phi_xy.*lnn_y.*Pi_y + phi_x.*lnn_yy.*Pi_y + phi_x.*lnn_y.*Pi_yy)
    _j = -ad*(phi_yy.*lnn_x.*Pi_y + phi_y.*lnn_xy.*Pi_y + phi_y.*lnn_x.*Pi_yy)
    _k = -eg*G_c - Pe_c - Pi_c + j_s

    return +(_a, _b, _c, _d, _e, _f, _g, _h, _i, _j, _k)
end


function A_partial_t(phi_x, phi_y, phi_s, j, jn_x, jn_y, jn_s, ue, n, Te, de2, ev, ad, am, eta)

    _a = de2 * (phi_x .* jn_y - phi_y .* jn_x)
    _b = de2 * ue.*jn_s / ev
    _c = (phi_s - ad*Pe_s ./ n) / am
    _d = eta*(j ./ Te.^(3/2))

    return +(_a, _b, _c, _d)
end


function compute_j_bohm(Te, Ti, n, phi, ev, ad)

    cs = sqrt.(Te + Ti)
    lambda = 2.695
    return ev * n.*cs.*(1 .- exp.(lambda .- phi ./(ad*Te)))
end


function compute_G(Ti, phi_c, Pi_c, u_s, er, ev, ad)

    return Ti.^(5/2) .* (er*(phi_c + ad*Pi_c) - 4*ev*u_s)
end


function compute_ue(u, jn, ev)

    return u - jn / ev
end


function density_source()

    return 0 #XXX no sources for now
end


function heat_source()

    return 0 #XXX no sources for now
end


function vorticity_lhs(n, lnn_x, lnn_y, L, Dx, Dy, P0, P1, P2, P3, R1, R2, R3)

    N = Diagonal(n)
    lnN_x = Diagonal(lnn_x)
    lnN_y = Diagonal(lnn_y)

    _A = P0 * N * (L + lnN_x*Dx + lnN_y*Dy)
    _B = P1  # XXX try this condition without the reflection at first
    _C = 0.5*P2*(I + R2)
    _D = 0.5*P3*(I + R3)

    return +(_A, _B, _C, _D)
end


function vorticity_rhs()

    lambda = 2.695
    _a = P0*(w - ad*(Pi_xx + Pi_yy))
    _b = P1*lnn_b # XXX try this condition without the reflection at first
    _c = 0.5*lambda*P2*(I + R2)*lnTe
    _d = 0.5*lambda*P3*(I + R3)*lnTe

    return +(_a, _b, _c, _d)
end


function diffusion_neumann_lhs()

    _A = P0 * (I - k*L)
    _B = P1*(I - R1)
    _C = P2*(I - R2)
    _D = P3*(I - R3)

    return +(_A, _B, _C, _D)
end


function diffusion_dirichlet_lhs()

    _A = P0 * (I - k*L)
    _B = 0.5*P1*(I + R1)
    _C = 0.5*P2*(I + R2)
    _D = 0.5*P3*(I + R3)

    return +(_A, _B, _C, _D)
end


function diffusion_u_lhs()

    _A = P0 * (I - k*L)
    _B = P1*(I - R1)
    _C = P2*(I - R2)
    _D = 0.5*P3*(I + R3)

    return +(_A, _B, _C, _D)
end


function diffusion_u_rhs()

    return P0*u + 0.5*P3*(I + R3)*sqrt.(Te + Ti)
end


function diffusion_A_rhs()

    return P0*A + 0.5*P3*(I + R3)*Abohm
end


function penalized_forward_euler(f, f_t, dt, p0, p1, p2, p3, q1, q2, q3)

    c1 = p1 / q1
    c2 = p2 / q2
    c3 = p3 / q3

    _a = 1 .+ dt*(c1 + c2 + c3)
    _b = +(p0.*f_t, c1.*f1, c2.*f2, c3.*f3)

    return (f + dt*_b) ./ _a
end





























#
