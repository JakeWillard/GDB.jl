

taking_derivatives_in_serial = quote

    phi_x = Dx*phi[:,2]
    phi_y = Dy*phi[:,2]
    phi_xx = Dxx*phi[:,2]
    phi_yy = Dyy*phi[:,2]
    phi_xy = Dxy*phi[:,2]

    psi_x = Dx*psi[:,2]
    psi_y = Dy*psi[:,2]
    psi_xx = Dxx*psi[:,2]
    psi_yy = Dyy*psi[:,2]
    psi_xy = Dxy*psi[:,2]

    lnn_x = Dx*lnn[:,2]
    lnn_y = Dy*lnn[:,2]
    lnn_xx = Dxx*lnn[:,2]
    lnn_yy = Dyy*lnn[:,2]
    lnn_xy = Dxy*lnn[:,2]

    lnTe_x = Dx*lnTe
    lnTe_y = Dy*lnTe
    lnTe_xx = Dxx*lnTe
    lnTe_yy = Dyy*lnTe
    lnTe_xy = Dxy*lnTe

    lnTi_x = Dx*lnTi
    lnTi_y = Dy*lnTi
    lnTi_xx = Dxx*lnTi
    lnTi_yy = Dyy*lnTi
    lnTi_xy = Dxy*lnTi

    w_x = Dx*w
    w_y = Dy*w
    w_xx = Dxx*w
    w_yy = Dyy*w
    w_xy = Dxy*w

    u_x = Dx*u
    u_y = Dy*u
    u_xx = Dxx*u
    u_yy = Dyy*u
    u_xy = Dxy*u

    A_x = Dx*A
    A_y = Dy*A
    A_xx = Dxx*A
    A_yy = Dyy*A
    A_xy = Dxy*A

end
