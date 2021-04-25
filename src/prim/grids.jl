

struct Grid

    points::Matrix{Float64}
    projection::SparseMatrixCSC
    inverse_projection::SparseMatrixCSC
    Nx::Int32
    Ny::Int32
    Nk::Int32
    x_stencil::Matrix{Float64}
    y_stencil::Matrix{Float64}
    dx::Float64
    dy::Float64

end


function stencil(m::Int)

    ic = Int(ceil(m/2.0))
    M = zeros((m, m))
    for i=1:m
        for j=1:m
            M[i, j] = (i - ic)^(j-1) / factorial(j-1)
        end
    end
    return inv(M)
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

    # compute stencils
    x_stencil = stencil(mx)
    y_stencil = stencil(my)


    return Grid(points, P, Pinv, Nx, Ny, Nk, x_stencil, y_stencil, dx, dy)
end
