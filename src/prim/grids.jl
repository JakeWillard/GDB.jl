

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


function vec_to_mesh(vec::Vector{Float64}, grid::Grid)

    vals = grid.inverse_projection * vec

    Nx = grid.Nx
    Ny = grid.Ny
    out = zeros(Float64, (Nx, Ny))

    for j=1:Ny
        out[:,j] = vals[1+(j-1)*Nx:j*Nx]
    end

    return out
end


function f_to_grid(f::Function, grid::Grid)

    Nk = grid.Nk
    vec = zeros(Float64, Nk)

    for k=1:Nk
        x = grid.points[1,k]
        y = grid.points[2,k]
        vec[k] = f(x, y)
    end

    return vec
end

function load_matrix(group)

    V = group["V"][:]
    Is = group["I"][:]
    J = group["J"][:]
    n, m = group["size"][:]
    return sparse(Is, J, V, n, m)
end

function save_matrix(group, A::SparseMatrixCSC)

    n, m = size(A)
    Is, J, V = findnz(A)
    group["I"] = Int32[Is...]
    group["J"] = Int32[J...]
    group["V"] = Float64[V...]
    group["size"] = Int[n, m]
end

function save_grid(fid, grd::Grid)

    create_group(fid, "Grid/Proj")
    create_group(fid, "Grid/ProjInv")
    save_matrix(fid["Grid/Proj"], grd.projection)
    save_matrix(fid["Grid/ProjInv"], grd.inverse_projection)

    fid["Grid/points"] = grd.points[:,:]
    fid["Grid/dims"] = Float64[grd.Nx, grd.Ny, grd.Nk, grd.mx, grd.my]
    fid["Grid/spacing"] = Float64[grd.dx, grd.dy]
    fid["Grid/Mxinv"] = grd.Mxinv[:,:]
    fid["Grid/Myinv"] = grd.Myinv[:,:]
    fid["Grid/MxyinvT"] = grd.MxyinvT[:,:]

end


function load_grid(fid)

    P = load_matrix(fid["Grid/Proj"])
    Pinv = load_matrix(fid["Grid/ProjInv"])
    points = fid["Grid/points"][:,:]
    Nx, Ny, Nk, mx, my = fid["Grid/dims"][:]
    dx, dy = fid["Grid/spacing"][:]
    Mxinv = fid["Grid/Mxinv"][:,:]
    Myinv = fid["Grid/Myinv"][:,:]
    MxyinvT = fid["Grid/MxyinvT"][:,:]

    return Grid(points, P, Pinv, mx, my, Nx, Ny, Nk, Mxinv, Myinv, MxyinvT, dx, dy)
end
