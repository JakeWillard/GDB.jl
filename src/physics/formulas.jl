
plasma_frequency(n, q, m) = sqrt(4*pi*n*q^2 / m)
gyro_frequency(q, B, m) = q*B / (m * 3e10)
electron_collision_time(n, T, q, m, l) = 3*sqrt(m)*T^(3/2) / (4*sqrt(2*pi)*l*n*q^4)
ion_collision_time(n, T, q, m, l) = 3*sqrt(m)*T^(3/2) / (4*sqrt(pi)*l*n*q^4)
sound_speed(T, m) = sqrt(T / m)
plasma_beta(n, T, B) = 4*pi*n*T / B^2


#NOTE: input cgs units, except T0, input eV for T0
function dimensionless_parameters(a, R0, n0, T0, B0)

    # constants
    c = 3.00e10
    q = 4.8032e-10
    me = 9.1094e-28
    mp = 1.6726219e-24

    # convert T0 from eV to erg
    T0 = T0 * 1.60218e-12

    # assume deuterium ions
    mi = 2 * mp

    # important parameters
    wpe = plasma_frequency(n0, q, me)
    wci = gyro_frequency(q, B0, mi)
    tau_i = ion_collision_time(n0, T0, q, mi, 1)
    tau_e = electron_collision_time(n0, T0, q, me, 1)
    beta0 = plasma_beta(n0, T0, B0)
    cs0 = sound_speed(T0, mi)
    t0 = sqrt(a*R0/2) / cs0

    # parameters for gdb model
    am = R0*beta0/a
    ad = cs0^2 * t0 / (wci*a^2)
    ki = 2*3.9*t0*tau_i*T0 / (R0^2 * mi)
    ke = 2*3.2*t0*tau_e*T0 / (R0^2 * me)
    er = 2*a/R0
    eg = 0.08*tau_i / t0
    ev = cs0*t0/R0
    de2 = c / (a*wpe)
    eta = 0.51*t0*de2 / tau_e

    return Float64[am, ad, ki, ke, er, eg, ev, de2, eta]
end


function partial_s(f, f_x, f_y, psi_x, psi_y, Ds, am)

    return Ds*f + am*(psi_x .* f_y - psi_y .* f_x)
end


function partial_ss(f, f_s, f_x, f_y, f_xy, f_xx, f_yy, psi_x, psi_y, psi_xy, psi_xx, psi_yy, Dss, Ds, Dx, Dy, am)

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
    return k*T.^(5/2) .* _a ./ n
end


function u_total_ion_derivative(Pe_s, Pi_s, G_s, n, Ti, u_c, ev, eg, er, ad)

    _a = -ev*(Pe_s + Pi_s) ./ n
    _b = -4*ev*eg*G_s ./ n
    _c = er*ad*(Ti .* u_c)

    return +(_a, _b, _c)
end


function general_partial_derivative(f_Dt, f_x, f_y, f_s, phi_x, phi_y, v)

    return f_Dt - phi_x.*f_y + phi_y.*f_x - v.*f_s
end



function w_partial_t(n, lnn_x, lnn_y, lnn_xy, lnn_xx, lnn_yy,
                     Pe_x, Pe_y, Pe_xy, Pe_xx, Pe_yy, Pe_xxx, Pe_yyy, Pe_xxy, Pe_xyy, Pe_c,
                     phi_x, phi_y, phi_xy, phi_xx, phi_yy, phi_xxx, phi_yyy, phi_xxy, phi_xyy,
                     G_c, Pi_c, j_s, ad, eg)

    _a = (lnn_x.*lnn_y + lnn_xy).*(Pe_y.*phi_y - Pe_x.*phi_x)
    _b = (lnn_x.^2 + lnn_xx).*Pe_x.*phi_y - (lnn_y.^2 + lnn_yy).*Pe_y.*phi_x
    _c = -Pe_x.*lnn_x.^2 .*phi_y + Pe_x.*lnn_x.*lnn_y.*phi_x
    _d = Pe_x.*lnn_x.*phi_xy - Pe_x.*lnn_y.*phi_xx + Pe_xx.*lnn_x.*phi_y
    _e = -Pe_xx.*lnn_y.*phi_x - Pe_xx.*phi_xy - Pe_xxx.*phi_y
    _f = -Pe_y.*lnn_x.*lnn_y.*phi_y + Pe_y.*lnn_x.*phi_yy
    _g = Pe_y.*lnn_y.^2 .*phi_x - Pe_y.*lnn_y.*phi_xy + Pe_yy.*lnn_x.*phi_y
    _h = -Pe_yy.*lnn_y.*phi_x + Pe_yy.*phi_xy + Pe_yyy.*phi_x + phi_x.*Pe_xxy
    _i = phi_xx.*Pe_xy - ad.*phi_y.*Pe_xyy - phi_yy.*Pe_xy
    _j = lnn_x.*phi_x.*phi_xy - lnn_x.*phi_xx.*phi_y + lnn_y.*phi_x.*phi_yy - lnn_y.*phi_y.*phi_xy
    _k = -lnn_x.*phi_xx.*phi_y + lnn_y.*phi_x.*phi_yy - lnn_y.*phi_y.*phi_xy
    _l = phi_x.*phi_yyy + phi_x.*phi_xxy - phi_xxx.*phi_y - phi_y.*phi_xyy

    _A = ad * +(_a, _b, _c, _d, _e, _f, _g, _h, _i)
    _B = n .* +(_j, _k, _l)

    return _A + _B - eg*G_c - Pe_c - Pi_c + j_s
end


function A_partial_t(phi_x, phi_y, phi_s, j, jn, jn_x, jn_y, jn_s, Pe_s, u, n, Te, de2, ev, ad, am, eta)

    _a = de2 * (phi_x .* jn_y - phi_y .* jn_x)
    _b = de2 * (u - jn / ev).*jn_s / ev
    _c = (phi_s - ad*Pe_s ./ n) / am
    _d = eta*(j ./ Te.^(3/2))

    return +(_a, _b, _c, _d)
end


function Abohm_def(Te, Ti, n, phi, er, ev, ad, de2)

    cs = sqrt.(Te + Ti)
    lambda = 2.695
    jbohm = ev * n.*cs.*(1 .- exp.(lambda .- phi ./(ad*Te))) / (ad*er)

    return -de2*jbohm
end


function G_def(Ti, phi_c, Pi_c, u_s, er, ev, ad)

    return Ti.^(5/2) .* (er*(phi_c + ad*Pi_c) - 4*ev*u_s)
end
