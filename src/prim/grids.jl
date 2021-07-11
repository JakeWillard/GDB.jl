

struct Grid

    points::Matrix{Float64}
    origin::Vector{Float64}
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


function Grid(walls::Vector{Wall}, deltas::Vector{Float64}, corners::Matrix{Float64}, N::Vector{Int}, m::Vector{Int})

    Nx, Ny = N
    Nxy = Nx*Ny
    points = zeros(Float64, (2, Nxy))

    # c1 is bottom left corner of bounding box, c2 is top right
    c1 = corners[:,1]
    c2 = corners[:,2]
    dx = (c2[1] - c1[1]) / Nx
    dy = (c2[2] - c1[2]) / Ny
    Nk = 0

    proj_rows = Int32[i for i=1:Nxy]
    proj_cols = zeros(Int32, Nxy)
    proj_vals = ones(Int32, Nxy)

    for j=1:Ny
        for i=1:Nx
            x = c1[1] + (i-1) * dx
            y = c1[2] + (j-1) * dy
            inside = true

            for l=1:length(walls)
                if smoothstep(x, y, deltas[l], walls[l]) == 0
                    inside = false
                    break
                end
            end

            if inside
                Nk += 1
                points[:,Nk] = [x, y]
                proj_cols[Nk] = i + Nx*(j-1)
            end

        end
    end

    P = sparse(proj_rows[1:Nk], proj_cols[1:Nk], proj_vals[1:Nk], Nk, Nxy)
    Pinv = transpose(P)

    # compute stencils
    mx, my = m
    Mxinv = stencil1d(mx)
    Myinv = stencil1d(my)
    MxyinvT = stencil2d(mx, my)

    return Grid(points[:,1:Nk], c1, P, Pinv, mx, my, Nx, Ny, Nk, Mxinv, Myinv, MxyinvT, dx, dy)
end


function Grid(grd::Grid, walls::Vector{Wall}, deltas::Vector{Float64}, corners::Matrix{Float64}, n::Vector{Int}, m::Vector{Int})

    nx, ny = n
    Nxy = nx*ny*grd.Nk
    points = zeros(Float64, (2, Nxy))

    c1 = corners[:,1]
    c2 = corners[:,2]
    dx = grd.dx / nx
    dy = grd.dy / ny
    Nk = 0

    proj_rows = Int32[i for i=1:Nxy]
    proj_cols = zeros(Int32, Nxy)
    proj_vals = ones(Int32, Nxy)

    for k=1:grd.Nk
        for j=1:ny
            for i=1:nx
                x = grd.points[1,k] + (i-1) * dx
                y = grd.points[2,k] + (j-1) * dx
                ii = Int(floor((x - c1[1]) / dx))
                jj = Int(floor((y - c1[2]) / dy))
                inside = true

                for l=1:length(walls)
                    if smoothstep(x, y, deltas[l], walls[l]) == 0
                        inside = false
                        break
                    end
                end

                if inside
                    Nk += 1
                    points[:,Nk] = [x, y]
                    proj_cols[Nk] = ii + nx*grd.Nx*(jj-1)
                end

            end
        end
    end

    P = sparse(proj_rows[1:Nk], proj_cols[1:Nk], proj_vals[1:Nk], Nk, nx*ny*grd.Nx*grd.Ny)
    Pinv = transpose(P)

    # compute stencils
    mx, my = m
    Mxinv = stencil1d(mx)
    Myinv = stencil1d(my)
    MxyinvT = stencil2d(mx, my)

    return Grid(points[:,1:Nk], c1, P, Pinv, mx, my, nx*grd.Nx, ny*grd.Ny, Nk, Mxinv, Myinv, MxyinvT, dx, dy)
end


function Grid(walls::Array{Wall, 1}, deltas::Vector{Float64}, corners::Matrix{Float64}, N::Int, m::Int)

    return Grid(walls, deltas, corners, Int[N, N], Int[m, m])
end


function Grid(grd::Grid, walls::Array{Wall, 1}, deltas::Vector{Float64}, corners::Matrix{Float64}, n::Int, m::Int)

    return Grid(grd, walls, deltas, corners, Int[n, n], Int[m, m])
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
