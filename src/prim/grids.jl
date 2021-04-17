

struct Grid

    points::Matrix{Float64}
    projection::SparseMatrixCSC
    inverse_projection::SparseMatrixCSC

end


function Grid(inside::Function, Nx::Int, Ny::Int)

    N = Nx * Ny
    Nk = 0
    points = zeros(Float64, (2, N))

    proj_rows = Int32[i for i=1:N]
    proj_cols = zeros(Int32, N)
    proj_vals = ones(Int32, N)

    for i=1:Nx
        for j=1:Ny
            x = i / Nx
            y = j / Nx
            if inside([x, y])
                Nk += 1
                points[:,Nk] = [x, y]
                proj_cols[Nk] = i + Nx*(j-1)
            end
        end
    end

    P = sparse(proj_rows[1:Nk], proj_cols[1:Nk], proj_vals[1:Nk], Nk, N)
    Pinv = transpose(P)

    return Grid(points, P, Pinv)

end


# might be nice to construct 2^n x 2^n sized grids
function Grid(inside::Function, n::Int)

    N = 2^n
    Nx = Int(sqrt(N))
    return Grid(inside, Nx, Nx)
end
