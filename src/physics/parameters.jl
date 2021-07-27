

plasma_frequency(n, q, m) = sqrt(4*pi*n*q^2 / m)

gyro_frequency(q, B, m) = q*B / (m * 3e8)

electron_collision_time(n, T, q, m, l) = 3*sqrt(m)*T^(3/2) / (4*sqrt(2*pi)*l*n*q^4)

ion_collision_time(n, T, q, m, l) = 3*sqrt(m)*T^(3/2) / (4*sqrt(pi)*l*n*q^4)

sound_speed(T, m) = sqrt(T / m)

plasma_beta(n, T, B) = 4*pi*n*T / B^2


#NOTE: input cgs units, except T0, inpute eV for T0
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
    ad = cs0^2 * t0 / (omega_ci*a^2)
    ki = 2*3.9*t0*tau_i*Ti0 / (R0^2 * mi)
    ke = 2*3.2*t0*tau_e*Te0 / (R0^2 * me)
    er = 2*a/R0
    eg = 0.08*tau_i / t0
    ev = cs0*t0/R0
    de2 = c / (a*omega_pe0)
    eta = 0.51*t0*de2 / tau_e

    return Float64[am, ad, ki, ke, er, eg, ev, de2, eta]
end
