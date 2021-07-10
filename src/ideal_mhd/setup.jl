


function Rectangle(Lx, Ly)

    vert = zeros(Float64, (2, 5))
    vert[:,1] = Float64[-Lx/2, -Ly/2]
    vert[:,2] = Float64[Lx/2, -Ly/2]
    vert[:,3] = Float64[Lx/2, Ly/2]
    vert[:,4] = Float64[-Lx/2, Ly/2]
    vert[:,5] = Float64[-Lx/2, -Ly/2]

    return PolyWall(vert)
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

    outer_wall = Rectangle(Lx, Ly)
    x_flx = FluxWall(psi, psi(Lx/2 - h, 0), bx, by, ds/10.0)
    y_flx = FluxWall(psi, psi(0, Ly/2 - h), bx, by, ds/10.0)

    # define corners and deltas
    corners = Float64[-Lx/2-h Lx/2+h; -Ly/2-h Ly/2+h]
    crs_deltas = Float64[0.5, -0.5, 0.5]
    fine_deltas = Float64[0.01, -0.01, 0.01]

    # make grids
    crs_grd = Grid([outer_wall, x_flx, y_flx], crs_deltas, corners, N, m)
    fine_grd = Grid(crs_grd, [outer_wall, x_flx, y_flx], fine_deltas, corners, n, m)

    return SimulationSetup(fine_grd)
end
