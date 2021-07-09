


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


function SimulationSetup(Lx, Ly, h, ds, N, m)

    psi(x,y) = Psi(x, y, Lx, Ly)
    bx(x,y) = b(x,y,Lx,Ly)[1]
    by(x,y) = b(x,y,Lx,Ly)[2]

    # make coarse grid
    psi_a = psi(Lx/2.0 - h - 5*ds, 0)
    psi_b = psi(0, Ly/2.0 - h - 5*ds)
    dp_a = psi(Lx/2.0 - h - 4*ds, 0) - psi_a
    dp_b = psi(0, Ly/2.0 - h - 4*ds) - psi_b
    println(psi_a)
    println(psi_b)

    outer_wall = Rectangle(Lx+5*ds, Ly+5*ds, ds)
    pos_flx = FluxWall(psi, psi_a, bx, by, ds, dp_a)
    neg_flx = FluxWall(psi, psi_b, bx, by, ds, dp_b)

    corners = Float64[-Lx/2-h Lx/2+h; -Ly/2-h Ly/2+h]

    crs_grd = Grid([outer_wall, pos_flx, neg_flx], corners, N, m)

    return SimulationSetup(crs_grd)
end
