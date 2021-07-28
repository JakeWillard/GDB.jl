using Base: Int32
struct Diffresgrids
    R1::Matrix{Float64}
    R2::Matrix{Float64}
    R3::Matrix{Float64}
    I1::Matrix{Float64}
    I2::Matrix{Float64}
    I3::Matrix{Float64}
    Nx::Int32
    Ny::Int32
    Nx1::Int32
    Ny1::Int32
    Nx2::Int32
    Ny2::Int32
    Nx3::Int32
    Ny3::Int32


end

function multigrid_grids(crs_grd::Grid, walls::Vector{Wall}, deltas::Vector{Float64}, corners::Matrix{Float64}, m::Vector{Int})
    # fine_grid = Grid(crs_grd, Rectangle, thickness of step (one element vector), x/y position of bounding box (larger) and top right corner, 2, 3+)
    
    # defining the new grids by a change of 2 in resolution
    fine_grid = Grid(crs_grd, walls, deltas, corners, 2, m)
    finer_grid = Grid(fine_grid, walls, deltas, corners, 2, m)
    finest_grid = Grid(finer_grid, walls, deltas, corners, 2, m)

    # defining projections at each of the resolutions
    P1 = crs_grd.projection
    P1inv = crs_grd.inverse_projection

    P2 = fine_grid.projection
    P2inv = fine_grid.inverse_projection

    P3 = finer_grid.projection
    P3inv = finer_grid.inverse_projection

    P4 = finest_grid.projection
    P4inv = finest_grid.inverse_projection

    # obtaining restriction/interpolation matrices in Cartesian coordinates at all resolutions
    # need to resolve difference in resolution between given grids and the resolutions 
    # calculated in the transformation_matrices function
    Nx = crs_grd.Nx
    Ny = crs_grd.Ny
    R1, I1, Nxnew, Nynew = transformation_matrices(Nx, Ny)

    Nx1 = fine_grid.Nx
    Ny1 = fine_grid.Ny
    R2, I2, Nx1new, Ny1new = transformation_matrices(Nx1, Ny1)

    Nx2 = finer_grid.Nx
    Ny2 = finer_grid.Ny
    R3, I3, Nx2new, Ny2new = transformation_matrices(Nx2, Ny2)

    # calculating generalized restriction/interpolation matrices
    R1p = P2 * R1 * P1inv
    R2p = P3 * R2 * P2inv
    R3p = P4 * R3 * P3inv

    I1p = P1 * I1 * P2inv
    I2p = P2 * I2 * P3inv
    I3p = P3 * I3 * P4inv

    return diffresgrids(R1p, R2p, R3p, I1p, I2p, I3p, Nx, Ny, Nx1, Ny1, Nx2, Ny2, Nx3, Ny3)
end