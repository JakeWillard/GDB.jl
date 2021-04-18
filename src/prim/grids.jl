

struct Grid

    points::Matrix{Float64}
    projection::SparseMatrixCSC
    inverse_projection::SparseMatrixCSC
    Dx::SparseMatrixCSC
    Dxx::SparseMatrixCSC
    Dxy::SparseMatrixCSC
    Dy::SparseMatrixCSC
    Dyy::SparseMatrixCSC
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


function Grid(inside::Function, Nx::Int, Ny::Int, m::Int)

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

    # compute derivative matrices
    stens = stencil(m)
    Ix = sparse(I, Nx, Nx)
    Iy = sparse(I, Ny, Ny)
    Dx_1d = sparse(zeros(Float64, (N, N)))
    Dxx_1d = sparse(zeros(Float64, (N, N)))
    Dy_1d = sparse(zeros(Float64, (N, N)))
    Dyy_1d = sparse(zeros(Float64, (N, N)))

    for i=1:m
        ic = Int(ceil(m/2.0))
        k = i - ic
        vec_x = stens[2,i] * ones(Nx - abs(k))
        vec_xs = stens[3,i] * ones(Nx - abs(k))
        vec_y = stens[2,i] * ones(Ny - abs(k))
        vec_yy = stens[3,i] * ones(Ny - abs(k))

        Dx_1d += spdiagm(k => vec_x)
        Dxx_1d += spdiagm(k => vec_xx)
        Dy_1d += spdiagm(k => vec_y)
        Dyy_1d += spdiagm(k => vec_yy)

    end

    Dx = kron(Iy, Dx_1d)
    Dxx = kron(Iy, Dxx_1d)
    Dy = kron(Dy_1d, Ix)
    Dyy = kron(Dyy_1d, Ix)
    Dxy = Dx * Dy

    return Grid(points, P, Pinv, Dx, Dxx, Dxy, Dy, Dyy, dx, dy)
end



# might be nice to construct 2^n x 2^n sized grids
function Grid(inside::Function, n::Int)

    N = 2^n
    Nx = Int(sqrt(N))
    return Grid(inside, Nx, Nx)
end
