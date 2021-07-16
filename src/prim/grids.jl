

struct Grid

    points::Matrix{Float64}
    origin::Vector{Float64}
    projection::SparseMatrixCSC
    inverse_projection::SparseMatrixCSC
    restrictions::Vector{SparseMatrixCSC}
    interpolations::Vector{SparseMatrixCSC}
    nanmask::Matrix{Float64}
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

    nanmask = fill(NaN, (Nx, Ny))

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
                nanmask[i,j] = 1.0
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

    return Grid(points[:,1:Nk], c1, P, Pinv, nanmask, mx, my, Nx, Ny, Nk, Mxinv, Myinv, MxyinvT, dx, dy)
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

    # get nonzero rows from coarse projection matrix
    dummy, J, dummmy = findnz(grd.projection)

    proj_rows = Int32[i for i=1:Nxy]
    proj_cols = zeros(Int32, Nxy)
    proj_vals = ones(Int32, Nxy)

    # nanmask is useful for data visualization, represent points outside the
    # domain with NaN
    nanmask = fill(NaN, (nx*grd.Nx, ny*grd.Ny))

    for k=1:grd.Nk
        for j=1:ny
            for i=1:nx
                x = grd.points[1,k] + (i-1) * dx
                y = grd.points[2,k] + (j-1) * dx

                # compute (i,j) coordinates on coarse grid
                crs_i = rem(J[k], grd.Ny)
                crs_j = div(J[k], grd.Ny)

                # compute (i,j) coordinates on fine grid
                fn_i = (crs_i-1)*nx + i
                fn_j = (crs_j-1)*ny + j

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
                    proj_cols[Nk] = fn_i + nx*grd.Nx*(fn_j-1)
                    nanmask[fn_i, fn_j] = 1.0
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

    return Grid(points[:,1:Nk], c1, P, Pinv, nanmask, mx, my, nx*grd.Nx, ny*grd.Ny, Nk, Mxinv, Myinv, MxyinvT, dx, dy)
end


function Grid(walls::Array{Wall, 1}, deltas::Vector{Float64}, corners::Matrix{Float64}, N::Int, m::Int)

    return Grid(walls, deltas, corners, Int[N, N], Int[m, m])
end


function Grid(grd::Grid, walls::Array{Wall, 1}, deltas::Vector{Float64}, corners::Matrix{Float64}, n::Int, m::Int)

    return Grid(grd, walls, deltas, corners, Int[n, n], Int[m, m])
end



# XXX: multigrid stuff


function restriction1d(N)

    R = zeros(Float64, (N, 2*N))
    for i=1:N-1
        R[i, 2*(i-1)+1:2*(i-1)+3] = Float64[1, 2, 1]
    end
    R[N,N-1:N] = Float64[2, 1]

    return R / 4.0
end


function intergrid_transforms(crs_grd, fine_grd)

    Rxc = restriction1d(crs_grd.Nx)
    Ryc = restriction1d(crs_grd.Ny)
    Rc = kron(Ryc, Rxc)
    Ic = 4*transpose(Rc)

    R = crs_grd.projection * Rc * fine_grd.inverse_projection
    I = fine_grd.projection * Ic * crs_grd.inverse_projection

    return R, I
end


function Grid(Nl::Int, walls::Vector{Wall}, deltas::Vector{Float64}, corners::Matrix{Float64}, N::Vector{Int}, m::Vector{Int})

    Nx_crs, Ny_crs = N
    grd_0 = Grid(walls, deltas, corners, N, m)
    grd_1 = Grid(grd_0, walls, deltas, corners, 2, m)
    R0, I0 = intergrid_transforms(grd_0, grd_1)

    restrictions = SparseMatrixCSC[R0]
    interpolations = SparseMatrixCSC[I0]
    grds = Grid[grd_0, grd_1]

    for l=2:Nl

        grd_l = Grid(grds[l], walls, deltas, corners, 2, m)
        Rl, Il = intergrid_transforms(grds[l], grd_l)

    end

    


end











# XXX



function vec_to_mesh(vec, grid::Grid)


    vals = grid.inverse_projection * vec

    Nx = grid.Nx
    Ny = grid.Ny
    out = zeros(Float64, (Nx, Ny))

    for j=1:Ny
        out[:,j] = vals[1+(j-1)*Nx:j*Nx]
    end

    return out .* grid.nanmask
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
