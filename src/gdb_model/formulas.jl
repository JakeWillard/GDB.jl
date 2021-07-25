

function compute_params()

   c = 1 #XXX
   cs0 = sqrt(T0 / mi)
   t0 = sqrt(a*R0/2) / cs0

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


function partial_s(f, f_x, f_y, psi_x, psi_y, am, Ds)

    return Ds*f + am*(psi_x .* f_y - psi_y .* f_x)
end


function partial_ss()


end


function lnn_total_ion_derivative()

   _a = er*phi_c - er*ad*Pe_c ./ n
   _b = er*ad*j_s ./ n - ev*vp_s
   return _a + _b
end


function lnTe_total_electron_derivative_adiabatic()

   _a = (2/3)*(lnn_Dt + er*ad*j.*lnn_s./n)
   _b = -(5/3)*er*ad*(Te .* lnTe_c)

   return _a + _b
end


function lnTe_partial_thermal()

end


function lnTi_total_ion_derivative_adiabatic()

   _a = (2/3)*lnn_Dt
   _b = (5/3)*er*ad*(Ti .* lnTi_c)

   return _a + _b
end


function lnTi_partial_thermal()

end


function u_total_ion_derivative()

   _a = -ev*(Pe_s + Pi_s) ./ n
   _b = -4*ev*eg*G_s ./ n
   _c = er*ad*(Ti .* u_c)

   return +(_a, _b, _c)
end


function general_partial_derivative(f_Dt, f_x, f_y, f_s, phi_x, phi_y, v)

   return f_Dt - phi_x.*f_y + phi_y.*f_x - v*f_s
end


function w_partial_t()

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


function A_partial_t()

   # poisson bracket term
   _a = de2 * (phi_x .* jn_y - phi_y .* jn_x)
   _b = de2 * ue.*jn_s / ev
   _c = (phi_s - ad*Pe_s ./ n) / am
   _d = eta*(j ./ Te.^(3/2))

   return +(_a, _b, _c, _d)
end


function compute_j_bohm()


end


function compute_G()

   return Ti.^(5/2) .* (er*(phi_c + ad*Pi_c) - 4*ev*u_s)
end


function compute_ue()


end


function density_source()


end


function heat_source()


end


function vorticity_lhs()

   N = Diagonal(n)
   lnN_x = Diagonal(lnn_x)
   lnN_y = Diagonal(lnn_y)
   return N * (L + lnN_x*Dx + lnN_y*Dy)
end


function vorticity_rhs()

   return w - ad*(Pi_xx + Pi_yy)
end


function diffusion_neumann_lhs()

end


function diffusion_neumann_rhs()

end


function diffusion_dirichlet_lhs()

end


function diffusion_dirichlet_rhs()

end


function diffusion_A_lhs()

end


function diffusion_A_rhs()

end


function penalized_forward_euler(f, g, dt, c0, c1, c2, c3)

   lam = 1 .+ dt*(c1 + c2 + c3)
   g2 = +(c0.*g, c1.*f1, c2.*f2, c3.*f3)
   return (f + dt*g2) ./ lam
end





























#
