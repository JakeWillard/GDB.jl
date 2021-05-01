

struct Grid

    points::Matrix{Float64}
    projection::SparseMatrixCSC
    inverse_projection::SparseMatrixCSC
    mx::Int32
    my::Int32
    Nx::Int32
    Ny::Int32
    Nk::Int32
    Mxinv::Matrix{Float64}
    Myinv::Matrix{Float64}
    MxyinvT::Matrix{Float64}
    dx::Float64
    dy::Float64

end


function Grid(inside::Function, Nx::Int, Ny::Int, mx::Int, my::Int)

    N = Nx * Ny
    Nk = 0
    points = zeros(Float64, (2, N))
    dx = 1 / Nx
    dy = 1 / Ny

    # compute projection matrices
    proj_rows = Int32[i for i=1:N]
    proj_cols = zeros(Int32, N)
    proj_vals = ones(Int32, N)

    for j=1:Ny
        for i=1:Nx
            x = (i-1) / Nx
            y = (j-1) / Ny
            if inside(x, y)
                Nk += 1
                points[:,Nk] = [x, y]
                proj_cols[Nk] = i + Nx*(j-1)
            end
        end
    end

    P = sparse(proj_rows[1:Nk], proj_cols[1:Nk], proj_vals[1:Nk], Nk, N)
    Pinv = transpose(P)

    # compute stencils
    Mxinv = stencil1d(mx)
    Myinv = stencil1d(my)
    MxyinvT = stencil2d(mx, my)


    return Grid(points[:,1:Nk], P, Pinv, mx, my, Nx, Ny, Nk, Mxinv, Myinv, MxyinvT, dx, dy)
end
