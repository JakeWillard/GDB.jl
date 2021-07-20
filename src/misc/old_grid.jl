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
                crs_i = rem(J[k], grd.Nx)
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