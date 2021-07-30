

function partial_s(f, f_x, f_y, psi_x, psi_y, Ds, am)

    return Ds*f + am*(psi_x .* f_y - psi_y .* f_x)
end


function partial_ss(f, f_x, f_y, f_xy, f_xx, f_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, Dss, Ds, Dx, Dy, am)

    _a = Dss*f + am*Ds*(psi_x.*f_y - psi_y.*f_x) + am*(psi_x.*Dy*f_s - psi_y.*Dx*f_s)
    _b = (f_x.*psi_xy + f_xx.*psi_y - f_y.*psi_xx - psi_x.*f_xy).*psi_y
    _c = (f_x.*psi_yy - f_y.*psi_xy - f_yy.*psi_x + psi_y.*f_xy).*psi_x

    return _a + am^2 * (_b - _c)
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


n: 0
lnn: x, y, xy, xx, yy,
Pe: x, y, xy, xx, yy, xxx, yyy, xxy, xyy, c
phi: x, y, xy, xx, yy, yyy, xxy, xxx, xyy
G_c, Pi_c, j_s, ad, eg



function w_partial_t(n, lnn_x, lnn_y, lnn_xy, lnn_xx, lnn_yy,
                     Pe_x, Pe_y, Pe_xy, Pe_xx, Pe_yy, Pe_xxx, Pe_yyy, Pe_xxy, Pe_xyy, Pe_c,
                     phi_x, phi_y, phi_xy, phi_xx, phi_yy, phi_xxx, phi_yyy, phi_xxy, phi_xyy,
                     G_c, Pi_c, j_s, ad, eg)

    _a = (lnn_x.*lnn_y + lnn_xy).*(Pe_y.*phi_y - Pe_x.*phi_x)
    _b = (lnn_x.^2 + lnn_xx).*Pe_x.*phi_y - (lnn_y.^2 + lnn_yy)*Pe_y.*phi_x
    _c = -Pe_x.*lnn_x.^2 .*phi_y + Pe_x.*lnn_x.*lnn_y.*phi_x
    _d = Pe_x.*lnn_x.*phi_xy - Pe_x.*lnn_y.*phi_xx + Pe_xx.*lnn_x.*phi_y
    _e = -Pe_xx.*lnn_y.*phi_x - Pe_xx.*phi_xy - Pe_xxx.*phi_y
    _f = -Pe_y.*lnn_x.*lnn_y.*phi_y + Pe_y.*lnn_x.*phi_yy
    _g = Pe_y.*lnn_y.^2 .*phi_x - Pe_y.*lnn_y.^phi_xy + Pe_yy.*lnn_x.*phi_y
    _h = -Pe_yy.*lnn_y.*phi_x + Pe_yy.*phi_xy + Pe_yyy.*phi_x + phi_x.*Pe_xxy
    _i = phi_xx.*Pe_xy - ad.*phi_y.*Pe_xyy - phi_yy.*Pe_xy
    _j = lnn_x.*phi_x*phi_xy - lnn_x.*phi_xx.*phi_y + lnn_y.*phi_x.*phi_yy - lnn_y.*phi_y.*phi_xy
    _k = -lnn_x.*phi_xx.*phi_y + lnn_y.*phi_x.*phi_yy - lnn_y.*phi_y.*phi_xy
    _l = phi_x.*phi_yyy + phi_x.*phi_xxy - phi_xxx.*phi_y - phi_y.*phi_xyy

    _A = ad * +(_a, _b, _c, _d, _e, _f, _g, _h, _i)
    _B = n .* +(_j, _k, _l)

    return _A + _B - eg*G_c - Pe_c - Pi_c + j_s
end


function A_partial_t(phi_x, phi_y, phi_s, j, jn_x, jn_y, jn_s, ue, n, Te, de2, ev, ad, am, eta)

    _a = de2 * (phi_x .* jn_y - phi_y .* jn_x)
    _b = de2 * ue.*jn_s / ev
    _c = (phi_s - ad*Pe_s ./ n) / am
    _d = eta*(j ./ Te.^(3/2))

    return +(_a, _b, _c, _d)
end


function jbohm_def(Te, Ti, n, phi, ev, ad)

    cs = sqrt.(Te + Ti)
    lambda = 2.695
    return ev * n.*cs.*(1 .- exp.(lambda .- phi ./(ad*Te)))
end


function G_def(Ti, phi_c, Pi_c, u_s, er, ev, ad)

    return Ti.^(5/2) .* (er*(phi_c + ad*Pi_c) - 4*ev*u_s)
end


function ue_def(u, jn, ev)

    return u - jn / ev
end


# function density_source()
#
#     return 0 #XXX no sources for now
# end
#
#
# function heat_source()
#
#     return 0 #XXX no sources for now
# end


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
