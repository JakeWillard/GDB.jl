


function Rectangle(Lx, Ly, ds)

    vert = zeros(Float64, (2, 5))
    vert[:,5] = Float64[-Lx/2, -Ly/2]
    vert[:,4] = Float64[Lx/2, -Ly/2]
    vert[:,3] = Float64[Lx/2, Ly/2]
    vert[:,2] = Float64[-Lx/2, Ly/2]
    vert[:,1] = Float64[-Lx/2, -Ly/2]

    return PolyWall(vert, ds)
end


function Psi(x, y, Lx, Ly)

    return (2*x/Lx)^2 - (2*y/Ly)^2
end


function b(x, y, Lx, Ly)

    Bx = 4*y / Ly^2
    By = 4*x / Lx^2
    B = Float64[Bx, By]
    return B / norm(B)
end


struct SimulationSetup

    grd::Grid

end


function SimulationSetup(Lx, Ly, h, ds, N, n, m)

    psi(x,y) = Psi(x, y, Lx, Ly)
    bx(x,y) = b(x,y,Lx,Ly)[1]
    by(x,y) = b(x,y,Lx,Ly)[2]

    # define dimensions for coarse grid
    Lxc, Lyc = [Lx, Ly] .+ 5*ds
    psi_ac = psi(Lxc/2 - h, 0)
    psi_bc = psi(0, Lyc/2 - h)
    dp_ac = psi(Lxc/2 - h + ds, 0) - psi_ac
    dp_bc = psi(0, Lyc/2 - h + ds) - psi_bc

    # make walls for coarse grid
    outer_wall_c = Rectangle(Lxc, Lyc, ds)
    x_flx_c = FluxWall(psi, psi_ac, bx, by, ds, dp_ac)
    y_flx_c = FluxWall(psi, psi_bc, bx, by, ds, dp_bc)

    # make coarse grid
    corners = Float64[-Lxc/2-h Lxc/2+h; -Lyc/2-h Lyc/2+h]
    crs_grd = Grid([outer_wall_c, x_flx_c, y_flx_c], corners, N, m)

    # define dimensions for fine grid
    psi_af = psi(Lx/2 - h, 0)
    psi_bf = psi(0, Ly/2 - h)
    dp_af = psi(Lx/2 - h + ds, 0) - psi_af
    dp_bf = psi(0, Ly/2 - h + ds) - psi_bf

    # make fine grid
    outer_wall_f = Rectangle(Lx, Ly, ds)
    x_flx_f = FluxWall(psi, psi_af, bx, by, ds, dp_af)
    y_flx_f = FluxWall(psi, psi_bf, bx, by, ds, dp_bf)
    fine_grd = Grid(crs_grd, [outer_wall_f, x_flx_f, y_flx_f], corners, n, m)

    return SimulationSetup(fine_grd)
end
