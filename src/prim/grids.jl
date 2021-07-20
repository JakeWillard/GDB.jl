
function intergrid_1d_interp(N, m)

    ic = Int(ceil(m/2.0))
    Minv = stencil1d(m)
    sten = transpose(Minv) * Float64[(0.5)^(j-1)/factorial(j-1) for j=1:m]

    Ih = zeros(Float64, (N,N))
    Is = I + Ih
    for i=1:m
        k = i - ic
        vec = sten[i] * ones(N - abs(k))
        Ih += diagm(k => vec)
    end

    Iout = zeros(Float64, (2*N,N))
    Iout[1:2:end,:] = Is
    Iout[2:2:end,:] = Ih
    return Iout
end


function intergrid_transforms(Nx0, Ny0, mx, my, P0, P1)

    Ix = intergrid_1d_interp(Nx0, mx)
    Iy = intergrid_1d_interp(Ny0, my)
    Ic = kron(sparse(Iy), sparse(Ix))
    Rc = transpose(Ic) / 4.0

    Restr = P0 * Rc * transpose(P1)
    Interp = P1 * Ic * transpose(P0)
    return Restr, Interp
end


function edge_trim_box(walls::Vector{Wall}, deltas::Vector{Float64}, Nx::Int64, Ny::Int64, bottom_left::Vector{Float64}, top_right::Vector{Float64})

    deltaX = top_right[1] - bottom_left[1]
    deltaY = top_right[2] - bottom_left[2]
    dx = deltaX / Nx
    dy = deltaY / Ny

    points = zeros(Float64, (2, Nx*Ny))
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
    return points[:,1:Nk], P, nanmask, dx, dy
end


function edge_trim_double(points0::Matrix{Float64}, walls::Vector{Wall}, deltas::Vector{Float64}, dx0::Float64, dy0::Float64, Nx0::Int64, Ny0::Int64, bottom_left::Vector{Float64})

    dx = dx0 / 2.0
    dy = dy0 / 2.0
    Nk0 = size(points0)[2]

    points = zeros(Float64, (2, 4*Nk0))
    proj_rows = Int32[i for i=1:4*Nk0]
    proj_cols = zeros(Int32, 4*Nk0)
    proj_vals = ones(Int32, 4*Nk0)
    nanmask = fill(NaN, (2*Nx0, 2*Ny0))
    Nk = 0

    for l=1:Nk0
        for j=1:2
            for i=1:2
                x = points0[1,l] + (i-1)*dx
                y = points0[2,l] + (j-1)*dy
                deltaX = x - bottom_left[1]
                deltaY = y - bottom_left[2]
                ic = Int(round(deltaX/dx))
                jc = Int(round(deltaY/dy))

                if inside(x, y, deltas, walls)
                    Nk += 1
                    points[:,Nk] = [x, y]
                    proj_cols[Nk] = ic + 2*Nx0*(jc-1)
                    nanmask[ic, jc] = 1.0
                end

            end
        end
    end

    P = sparse(proj_rows[1:Nk], proj_cols[1:Nk], proj_vals[1:Nk], Nk, 4*Nx0*Ny0)
    return points[:,1:Nk], P, nanmask, dx, dy
end



#
#
#
#
#
#
# function edge_trim_double(Nx0::Int64, Ny0::Int64, dx0::Float64, dy0::Float64, points0::Matrix{Float64}, proj_cols0::Vector{Int32}, deltas::Vector{Float64}, walls::Vector{Wall})
#
#     dx = dx0 / 2.0
#     dy = dy0 / 2.0
#     Nk0 = size(points0)[2]
#
#     points = zeros(Float64, (2, 4*Nk0))
#     proj_rows = Int32[i for i=1:4*Nk0]
#     proj_cols = zeros(Int32, 4*Nk0)
#     proj_vals = ones(Int32, 4*Nk0)
#     nanmask = fill(NaN, (2*Nx0, 2*Ny0))
#     Nk = 0
#
#     for l=1:Nk0
#         for j=1:2
#             for i=1:2
#                 x = points0[1,l] + (i-1)*dx
#                 y = points0[2,l] + (j-1)*dy
#
#                 i0 = rem(proj_cols0[l], Nx0)
#                 j0 = div(proj_cols0[l], Ny0)
#                 i2 = (i0-1)*2 + i
#                 j2 = (j0-1)*2 + j
#
#                 if inside(x, y, deltas, walls)
#                     Nk += 1
#                     points[:,Nk] = [x, y]
#                     proj_cols[Nk] = i2 + 2*Nx0*(j2-1)
#                     nanmask[i2, j2] = 1.0
#                 end
#
#             end
#         end
#     end
#
#     P = sparse(proj_rows[1:Nk], proj_cols[1:Nk], proj_vals[1:Nk], Nk, 4*Nx0*Ny0)
#     return points[:,1:Nk], P, nanmask, dx, dy
# end



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

function Grid(walls::Vector{Wall}, deltas::Vector{Float64}, corners::Matrix{Float64}, N_levels::Int64, N0::Vector{Int}, m::Vector{Int})

    Nx, Ny = N0
    mx, my = m

    # trim the bounding box, then double the resolutions and compute Restr and Interp
    points, P0, nanmask, dx, dy = edge_trim_box(walls, deltas, Nx, Ny, corners[:,1], corners[:,2])
    points, P1, nanmask, dx, dy = edge_trim_double(points, walls, deltas/2, dx, dy, Nx, Ny, corners[:,1])
    Restr, Interp = intergrid_transforms(Nx, Ny, mx+1, my+1, P0, P1)
    Nx = Nx * 2
    Ny = Ny * 2

    restrictions = SparseMatrixCSC[Restr]
    interpolations = SparseMatrixCSC[Interp]

    for l=3:N_levels

        P0 = P1
        points, P1, nanmask, dx, dy = edge_trim_double(points, walls, deltas/l, dx, dy, Nx, Ny, corners[:,1])
        Restr, Interp = intergrid_transforms(Nx, Ny, mx+1, my+1, P0, P1)
        prepend!(restrictions, [Restr])
        append!(interpolations, [Interp])
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


function Grid(walls::Vector{Wall}, deltas::Vector{Float64}, corners::Matrix{Float64}, N_levels::Int64, N0::Int64, m::Int64)

    return Grid(walls, deltas, corners, N_levels, [N0, N0], [m, m])
end


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
