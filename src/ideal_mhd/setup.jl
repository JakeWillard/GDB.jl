


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


struct Setup

    grd::Grid
    P1::SparseMatrixCSC
    P2::SparseMatrixCSC
    P3::SparseMatrixCSC
    R1::SparseMatrixCSC
    R2::SparseMatrixCSC
    R3::SparseMatrixCSC
    Dx::SparseMatrixCSC
    Dy::SparseMatrixCSC
    L::SparseMatrixCSC

end


function Setup(Lx, Ly, h, ds, N, n, m)

    psi(x,y) = Psi(x, y, Lx, Ly)
    bx(x,y) = b(x,y,Lx,Ly)[1]
    by(x,y) = b(x,y,Lx,Ly)[2]

    outer_wall = Rectangle(Lx, Ly)
    x_flx = FluxWall(psi, psi(Lx/2 - h, 0), bx, by, ds/10.0)
    y_flx = FluxWall(psi, psi(0, Ly/2 - h), bx, by, ds/10.0)

    # define corners and deltas
    corners = Float64[-Lx/2-h Lx/2+h; -Ly/2-h Ly/2+h]
    crs_deltas = Float64[1, -4*(Lx-h)/Lx^2, 4*(Ly-h)/Ly^2] * 0.5
    fine_deltas = Float64[1, -4*(Lx-h)/Lx^2, 4*(Ly-h)/Ly^2] * ds*2

    # make grids
    crs_grd = Grid([outer_wall, x_flx, y_flx], crs_deltas, corners, N, m)
    fine_grd = Grid(crs_grd, [outer_wall, x_flx, y_flx], fine_deltas, corners, n, m)

    # fine_grd = crs_grd

    # make penalization matrices
    P1 = penalization_matrix(-4*ds*(Lx-h)/Lx^2, x_flx, fine_grd)
    P2 = penalization_matrix(4*ds*(Ly-h)/Ly^2, y_flx, fine_grd)
    P3 = penalization_matrix(ds, outer_wall, fine_grd)

    # make reflection matrices
    R1 = reflection_matrix(-4*ds*(Lx-h)/Lx^2, x_flx, fine_grd)
    println("")
    println("computed R1")
    R2 = reflection_matrix(4*ds*(Ly-h)/Ly^2, y_flx, fine_grd)
    println("computed R2")
    R3 = reflection_matrix(ds, outer_wall, fine_grd)
    println("computed R3")

    # make operators
    Dx = x_derivative(1, fine_grd)
    Dy = y_derivative(1, fine_grd)
    L = laplacian(fine_grd)

    return Setup(fine_grd, P1, P2, P3, R1, R2, R3, Dx, Dy, L)
end
