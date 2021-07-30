

serial_calcs4formulas = quote

    phi_x = Dx*phi[:,t]
    phi_y = Dy*phi[:,t]
    phi_xx = Dxx*phi[:,t]
    phi_yy = Dyy*phi[:,t]
    phi_xy = Dxy*phi[:,t]

    psi_x = Dx*psi[:,t]
    psi_y = Dy*psi[:,t]
    psi_xx = Dxx*psi[:,t]
    psi_yy = Dyy*psi[:,t]
    psi_xy = Dxy*psi[:,t]

    lnn_x = Dx*lnn[:,t]
    lnn_y = Dy*lnn[:,t]
    lnn_xx = Dxx*lnn[:,t]
    lnn_yy = Dyy*lnn[:,t]
    lnn_xy = Dxy*lnn[:,t]

    lnTe_x = Dx*lnTe[:,t]
    lnTe_y = Dy*lnTe[:,t]
    lnTe_xx = Dxx*lnTe[:,t]
    lnTe_yy = Dyy*lnTe[:,t]
    lnTe_xy = Dxy*lnTe[:,t]

    lnTi_x = Dx*lnTi[:,t]
    lnTi_y = Dy*lnTi[:,t]
    lnTi_xx = Dxx*lnTi[:,t]
    lnTi_yy = Dyy*lnTi[:,t]
    lnTi_xy = Dxy*lnTi[:,t]

    w_x = Dx*w[:,t]
    w_y = Dy*w[:,t]
    w_xx = Dxx*w[:,t]
    w_yy = Dyy*w[:,t]
    w_xy = Dxy*w[:,t]

    u_x = Dx*u[:,t]
    u_y = Dy*u[:,t]
    u_xx = Dxx*u[:,t]
    u_yy = Dyy*u[:,t]
    u_xy = Dxy*u[:,t]

    A_x = Dx*A[:,t]
    A_y = Dy*A[:,t]
    A_xx = Dxx*A[:,t]
    A_yy = Dyy*A[:,t]
    A_xy = Dxy*A[:,t]

end
