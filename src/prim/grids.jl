

function intergrid_transforms(Nx0, Ny0, P0, P1)

    Rx = zeros(Float64, (Nx0, 2*Nx0))
    Ry = zeros(Float64, (Ny0, 2*Ny0))
    for i=1:Nx-1
        Rx[i, 2*(i-1)+1:2*(i-1)+3] = Float64[1, 2, 1] / 4.0
    end
    for i=1:Ny-1
        Ry[i, 2*(i-1)+1:2*(i-1)+3] = Float64[1, 2, 1] / 4.0

    Rc = kron(Ry, Rx)
    Ic = 4*transpose(Rc)

    Restr = P0 * Rc * transpose(P1)
    Interp = P1 * Ic * transpose(P0)
    return Restr, Interp
end


function edge_trim_box(walls::Vector{Wall}, deltas::Vector{Float64}, Nx::Int32, Ny::Int32, bottom_left::Vector{Float64}, bottom_right::Vector{Float64})

    deltaX = bottom_right[1] - bottom_left[1]
    deltaY = bottom_right[2] - bottom_left[2]
    dx = deltaX / Nx
    dy = deltaY / Ny

    points = zeroes(Float64, (2, Nx*Ny))
    proj_rows = Int32[i for i=1:Nx*Ny]
    proj_cols = zeros(Int32, Nx*Ny)
    proj_vals = ones(Int32, Nx*Ny)
    nanmask = fill(NaN, (Nx, Ny))
    Nk = 0

    for j=1:Ny
        for i=1:Nx

            x = bottom_left[1] + (i-1)*dx
            y = bottom_left[2] + (j-1)*dx
            if inside(x, y, deltas, walls)
                Nk += 1
                points[:,Nk] = [x, y]
                proj_cols[Nk] = i + Nx*(j-1)
                nanmask[i,j] = 1.0
            end
        end
    end

    P = sparse(proj_rows[1:Nk], proj_cols[1:Nk], proj_vals[1:Nk], Nk, Nx*Ny)
    return points, P, nanmask, dx, dy
end


function edge_trim_double(Nx0::Int32, Ny0::Int32, dx0::Float64, dy0::Float64, points0::Matrix{Float64}, proj_cols0::Vector{Float64}, deltas::Vector{Float64}, walls::Vector{Wall})

    dx = dx0 / 2.0
    dy = dy0 / 2.0
    Nk0 = size(points0)[2]

    points = zeroes(Float64, (2, 4*Nk0))
    proj_rows = Int32[i for i=1:4*Nk0]
    proj_cols = zeros(Int32, 4*Nk0)
    proj_vals = ones(Int32, 4*Nk0)
    nanmask = fill(NaN, (2*Nx, 2*Ny))
    Nk = 0

    for l=1:Nk0
        for j=1:2
            for i=1:2
                x = points0[1,k] + (i-1)*dx
                y = points0[2,k] + (j-1)*dy

                i0 = rem(proj_cols0[l], Nx0)
                j0 = div(proj_cols0[l], Ny0)
                i2 = (i0-1)*2 + i
                j2 = (j0-1)*2 + j

                if inside(x, y, deltas, walls)
                    Nk += 1
                    points[:,Nk] = [x, y]
                    proj_cols[Nk] = i2 + 2*Nx0*(j2-1)
                    nanmask[i2, j2] = 1.0
                end

            end
        end
    end

    P = sparse(proj_rows[1:Nk], proj_cols[1:Nk], proj_vals[1:Nk], Nk, 4*Nx*Ny)
    return points, P, nanmask, dx, dy
end



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


function Grid(walls::Vector{Wall}, deltas::Vector{Float64}, corners::Matrix{Float64}, N_levels::Int32, N0::Vector{Int}, m::Vector{Int})

    Nx, Ny = N
    mx, my = m

    # trim the bounding box, then double the resolutions and compute Restr and Interp
    points, P0, nanmask, dx, dy = edge_trim_box(walls, deltas, Nx, Ny, corners[:,1], corners[:,2])
    points, P1, nanmask, dx, dy = edge_trim_double(Nx, Ny, dx, dy, points, findnz(P0)[2], deltas, walls)
    Restr, Interp = intergrid_transforms(Nx, Ny, P0, P1)
    Nx = Nx * 2
    Ny = Ny * 2

    restrictions = SparseMatrixCSC[Restr]
    interpolations = SparseMatrixCSC[Interp]

    for l=3:N_levels

        P0 = P1
        points, P1, nanmask, dx, dy = edge_trim_double(Nx, Ny, dx, dy, points, findnz[P0][2], deltas, walls)
        Restr, Interp = intergrid_transforms(Nx, Ny, P0, P1)
        prepend!(restrictions, Restr)
        append!(interpolations, Interp)
        Nx = Nx * 2
        Ny = Ny * 2
    end

    origin = corners[:,1]
    projection = P1
    inverse_projection = transpose(P1)
    Nk = size(points)[2]
    Mxinv = stencil1d(mx)
    Myinv = stencil1d(my)
    MxyinvT = stencil2d(mx, my)

    return Grid(points, origin, projection, inverse_projection, restrictions, interpolations,
                nanmask, mx, my, Nx, Ny, Nk, Mxinv, Myinv, MxyinvT, dx, dy)
end


function Grid(walls::Vector{Wall}, deltas::Vector{Float64}, corners::Matrix{Float64}, N_levels::Int32, N0::Int32, m::Int32)

    return Grid(walls, deltas, corners, N_levels, [N0, N0], [m, m])
end


#
# function Grid(walls::Vector{Wall}, deltas::Vector{Float64}, corners::Matrix{Float64}, N::Vector{Int}, m::Vector{Int})
#
#     Nx, Ny = N
#     Nxy = Nx*Ny
#     points = zeros(Float64, (2, Nxy))
#
#     # c1 is bottom left corner of bounding box, c2 is top right
#     c1 = corners[:,1]
#     c2 = corners[:,2]
#     dx = (c2[1] - c1[1]) / Nx
#     dy = (c2[2] - c1[2]) / Ny
#     Nk = 0
#
#     proj_rows = Int32[i for i=1:Nxy]
#     proj_cols = zeros(Int32, Nxy)
#     proj_vals = ones(Int32, Nxy)
#
#     nanmask = fill(NaN, (Nx, Ny))
#
#     for j=1:Ny
#         for i=1:Nx
#             x = c1[1] + (i-1) * dx
#             y = c1[2] + (j-1) * dy
#             inside = true
#
#             for l=1:length(walls)
#                 if smoothstep(x, y, deltas[l], walls[l]) == 0
#                     inside = false
#                     break
#                 end
#             end
#
#             if inside
#                 Nk += 1
#                 points[:,Nk] = [x, y]
#                 proj_cols[Nk] = i + Nx*(j-1)
#                 nanmask[i,j] = 1.0
#             end
#
#         end
#     end
#
#     P = sparse(proj_rows[1:Nk], proj_cols[1:Nk], proj_vals[1:Nk], Nk, Nxy)
#     Pinv = transpose(P)
#
#     # compute stencils
#     mx, my = m
#     Mxinv = stencil1d(mx)
#     Myinv = stencil1d(my)
#     MxyinvT = stencil2d(mx, my)
#
#     return Grid(points[:,1:Nk], c1, P, Pinv, nanmask, mx, my, Nx, Ny, Nk, Mxinv, Myinv, MxyinvT, dx, dy)
# end
#
#
# function Grid(grd::Grid, walls::Vector{Wall}, deltas::Vector{Float64}, corners::Matrix{Float64}, n::Vector{Int}, m::Vector{Int})
#
#     nx, ny = n
#     Nxy = nx*ny*grd.Nk
#     points = zeros(Float64, (2, Nxy))
#
#     c1 = corners[:,1]
#     c2 = corners[:,2]
#     dx = grd.dx / nx
#     dy = grd.dy / ny
#     Nk = 0
#
#     # get nonzero rows from coarse projection matrix
#     dummy, J, dummmy = findnz(grd.projection)
#
#     proj_rows = Int32[i for i=1:Nxy]
#     proj_cols = zeros(Int32, Nxy)
#     proj_vals = ones(Int32, Nxy)
#
#     # nanmask is useful for data visualization, represent points outside the
#     # domain with NaN
#     nanmask = fill(NaN, (nx*grd.Nx, ny*grd.Ny))
#
#     for k=1:grd.Nk
#         for j=1:ny
#             for i=1:nx
#                 x = grd.points[1,k] + (i-1) * dx
#                 y = grd.points[2,k] + (j-1) * dx
#
#                 # compute (i,j) coordinates on coarse grid
#                 crs_i = rem(J[k], grd.Nx)
#                 crs_j = div(J[k], grd.Ny)
#
#                 # compute (i,j) coordinates on fine grid
#                 fn_i = (crs_i-1)*nx + i
#                 fn_j = (crs_j-1)*ny + j
#
#                 inside = true
#                 for l=1:length(walls)
#                     if smoothstep(x, y, deltas[l], walls[l]) == 0
#                         inside = false
#                         break
#                     end
#                 end
#
#                 if inside
#                     Nk += 1
#                     points[:,Nk] = [x, y]
#                     proj_cols[Nk] = fn_i + nx*grd.Nx*(fn_j-1)
#                     nanmask[fn_i, fn_j] = 1.0
#                 end
#
#             end
#         end
#     end
#
#     P = sparse(proj_rows[1:Nk], proj_cols[1:Nk], proj_vals[1:Nk], Nk, nx*ny*grd.Nx*grd.Ny)
#     Pinv = transpose(P)
#
#     # compute stencils
#     mx, my = m
#     Mxinv = stencil1d(mx)
#     Myinv = stencil1d(my)
#     MxyinvT = stencil2d(mx, my)
#
#     return Grid(points[:,1:Nk], c1, P, Pinv, nanmask, mx, my, nx*grd.Nx, ny*grd.Ny, Nk, Mxinv, Myinv, MxyinvT, dx, dy)
# end
#
#
# function Grid(walls::Array{Wall, 1}, deltas::Vector{Float64}, corners::Matrix{Float64}, N::Int, m::Int)
#
#     return Grid(walls, deltas, corners, Int[N, N], Int[m, m])
# end
#
#
# function Grid(grd::Grid, walls::Array{Wall, 1}, deltas::Vector{Float64}, corners::Matrix{Float64}, n::Int, m::Int)
#
#     return Grid(grd, walls, deltas, corners, Int[n, n], Int[m, m])
# end


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
