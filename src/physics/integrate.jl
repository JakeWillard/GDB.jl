
# NOTE: there aren't any functions in this file, just quoted expressions that would make
# these functions unreadable if they were pasted here.
include("./tedious.jl")


function subcycle_temperature()

end


function integrate(N, lnn, lnTe, lnTi, u, w, A, psi, phi, j, Sn, STe, STi, Dx, Dy, Dxx, Dyy, Dxy, params)

    Te = exp.(lnTe[:,2])
    Ti = exp.(Ti[:,2])
    n = exp.(lnn[:,2])
    jn = j[:,2] ./ n
    G = compute_G
    Pe = n .* Te
    Pi = n .* Ti

    # integrate N time steps
    for _=1:N

        # this expression is defined in the included file "./tedious.jl"
        t = 2
        eval(serial_calcs4formulas)



        # compute the total time derivatives
        lnn_Dt = lnn_total_ion_derivative(phi_c, Pe_c, n[:,2], j_s, u_s, er, ev, ad)
        lnTe_Dt = lnTe_total_electron_derivative_adiabatic(lnn_Dt, j[:,2], lnn_s, n[:,2], Te, lnTe_c, er, ad)
        lnTi_Dt = lnTi_total_ion_derivative_adiabatic(lnn_Dt, Ti, lnTi_c, er, ad)
        u_Dt = u_total_ion_derivative(Pe_s, Pi_s, n, G_s, n, Ti, u_c, ev, eg, er, ad)

        # compute partial time derivatives
        lnn_t = general_partial_derivative(lnn_Dt, lnn_x, lnn_y, lnn_s, phi_x, phi_y, u[:,2])
        lnTe_t = general_partial_derivative(lnTe_Dt, lnTe_x, lnTe_y, lnTe_s, phi_x, phi_y, ue)
        lnTi_t = general_partial_derivative(lnTi_Dt, lnTi_x, lnTi_y, lnTi_s, phi_x, phi_y, u[:,2])
        u_t = general_partial_derivative(u_Dt, u_x, u_y, u_s, phi_x, phi_y, u[:,2])
        w_t = 0 #XXX compute later
        A_t = A_partial_t(phi_x, phi_y, phi_s, j, jn_x, jn_y, jn_s, ue, n, Te, de2, ev, ad, am, eta)

        # compute predictor step
    end





end


# evaluate auxillary variables: n, Te, Ti, G, etc
# evaluate derivatives of things
# evaluate total time derivatives
# evaluate partial time derivatives
# compute predictor step
# solve linear problems
# evaluate auxillary again
# evaluate derivatives again
# evaluate total time derivatives again
# evaluate partial time derivatives again
# compute final time step
# solve linear problems


# evaluate auxillary variables:
