

function distance_to_segment(x, y, vert0::Vector{Float64}, vert1::Vector{Float64})

    xs, ys = [x, y] - vert0
    t = begin dv = vert1 - vert0; atan(dv[2], dv[1]) end

    u = xs * cos(t) + ys * sin(t)
    v = -xs * sin(t) + ys * cos(t)

    # the sign function returns 0 at v=0, which kind of messes things up, so have it default to 1 instead.
    sgn = sign(v)
    if sgn == 0
        sgn = 1
    end

    if u < 0
        return norm([x,y] - vert0) * sgn
    elseif 0 <= u < norm(vert1 - vert0)
        return v
    else
        return norm([x,y] - vert1) * sgn
    end
end


function distance_to_arm(x, y, vert0::Vector{Float64}, vert1::Vector{Float64}, vert2::Vector{Float64})

    # find shortest distance
    dist1 = distance_to_segment(x, y, vert0, vert1)
    dist2 = distance_to_segment(x, y, vert1, vert2)
    dist = minimum(abs.([dist1, dist2]))

    # get orientation of the bend
    r1 = vert1 - vert0
    r2 = vert2 - vert1
    bend = sign(r1[1]*r2[2] - r1[2]*r2[1])
    if bend == 0
        bend = 1
    end

    # get correct sign
    sgn = all(bend*[dist1, dist2] .> 0) ? bend : -bend

    return sgn * dist
end


struct Mirror

    verts :: Matrix{Float64}   # all of the vertices
    arms :: Array{Float64, 3}  # vertices in groups of 3

end
Base.getindex(M::Mirror, i) = Mirror(M.verts[:,i], M.arms[:,:,i])


function Mirror(chains::Vector{Matrix{Float64}})

    verts = hcat(chains...)
    arms_list = Array{Float64, 3}[]

    for chain in chains
        N = size(chain)[2]
        a = zeros(2, 3, N)
        a[:,:,1] = hcat(chain[:,end], chain[:,1:2])
        a[:,:,end] = hcat(chain[:,end-1:end], chain[:,1])
        for i=2:size(chain)[2]-1
            a[:,:,i] = chain[:,i-1:i+1]
        end
        append!(arms_list, [a])
    end

    arms = zeros(2, 3, size(verts)[2])
    arms[:,1,:] = hcat([a[:,1,:] for a in arms_list]...)
    arms[:,2,:] = hcat([a[:,2,:] for a in arms_list]...)
    arms[:,3,:] = hcat([a[:,3,:] for a in arms_list]...)

    return Mirror(verts, arms)
end


function distance_to_mirror(x, y, M::Mirror)

    # get the closest vertex, return distance to arm represented by that vertex
    _, k = findmin([norm([x,y] - M.verts[:,i]) for i=1:size(M.verts)[2]])
    vert0 = M.arms[:,1,k]
    vert1 = M.arms[:,2,k]
    vert2 = M.arms[:,3,k]

    return distance_to_arm(x, y, vert0, vert1, vert2), k
end


# function mirror_image(x, y, dx, dy, M::Mirror)
#
#     i = 0
#     j = 0
#     dist = 0
#     steps = 0
#
#     while distance_to_mirror(x, y, M)[1] < 0
#         hs = hcat([[distance_to_mirror(x+ii*dx, y+jj*dy, M)[1] for ii=-1:1] for jj=-1:1]...)
#         di, dj = Tuple(findmax(hs)[2]) .- 2
#         x += di*dx
#         y += dj*dy
#         i += di
#         j += dj
#         dist += sqrt(di^2 + dj^2)
#         steps += di^2 + dj^2
#     end
#
#     while (dist > 0) && (steps > 0)
#         hs = hcat([[distance_to_mirror(x+ii*dx, y+jj*dy, M)[1] for ii=-1:1] for jj=-1:1]...)
#         di, dj = Tuple(findmax(hs)[2]) .- 2
#         x += di*dx
#         y += dj*dy
#         i += di
#         j += dj
#         dist -= sqrt(di^2 + dj^2)
#         steps -= 1
#     end
#
#     return x, y, i, j
# end


function mirror_image(x0, y0, dx, dy, M::Mirror)

    H(r) = distance_to_mirror(r[1], r[2], M)[1]
    gradH(r) = ForwardDiff.gradient(H, r)

    x = x0
    y = y0
    i = 0.0
    j = 0.0
    dist = 0
    steps = 0

    while H([x,y]) < 0
        v = gradH([x,y])
        v = v / norm(v)
        @assert !any(isinf.(v))
        i += v[1]
        j += v[2]
        x = x0 + dx*round(i)
        y = y0 + dy*round(j)
        # step =
        steps += 1
        dist += norm(round.(v))
        # dist += sqrt(dstep)
    end

    steps = 2 * steps
    while (dist > 0) && (steps > 0)
        v = gradH([x,y])
        v = v / norm(v)
        i += v[1]
        j += v[2]
        x = x0 + dx*round(i)
        y = y0 + dy*round(j)
        # dstep = i^2 + j^2
        steps -= 1
        dist -= norm(round.(v))
        # dist -= sqrt(dstep)
    end

    return x, y, round(i), round(j)
end


function smoothstep(x, y, s0, ds, M::Mirror)

    u = (distance_to_mirror(x, y, M)[1] - s0) / ds
    return (u < 0) ? 0.0 : ((0 < u < 1) ? 3*u^2 - 2*u^3 : 1.0)
end







# GOES IN DIFFERENT FILE
struct GhostData

    Proj :: SparseMatrixCSC
    R :: SparseMatrixCSC
    flip_factors :: Matrix{Int64}

end


function GhostData(M::Mirror, grd::Grid)

    chnl = RemoteChannel(()->Channel{Bool}(), 1)
    p = Progress(grd.Nk, desc="Resolving mirrors... ")
    @async while take!(chnl)
        next!(p)
    end

    map_out = map(1:grd.Nk) do k

        _, proj_js, _ = findnz(grd.Proj)
        k_cart = proj_js[k]
        h, sindex = distance_to_mirror(grd.points[:,k]..., M)
        flip_row = ones(Int64, size(M.verts)[2])

        if h < 0

            flip_row[sindex] = -1

            j = Int(ceil(k_cart / grd._Nx))
            i = k_cart - (j - 1) * grd._Nx
            x, y, di, dj = mirror_image(grd.points[:,k]..., grd.dx, grd.dy, M)
            i += di
            j += dj
            k_cart = i + (j-1)*grd._Nx
            k = 0
        end

        put!(chnl, true)
        Int64[k, k_cart, flip_row...]
    end

    C = hcat(map_out...)

    proj_j = C[1,:]
    proj_j = proj_j[proj_j .!= 0]
    r_j = C[2,:]
    flip_factors = C[3:end,:]

    Mirror = sparse([1:grd.Nk...], r_j, ones(grd.Nk), grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
    Proj = sparse([1:length(proj_j)...], proj_j, ones(length(proj_j)), length(proj_j), grd.Nk)

    return GhostData(Proj, Mirror, flip_factors)
end


function flip_segments(gd::GhostData, bdry)

    R = gd.R
    fac = gd.flip_factors[bdry]
    for j=1:size(factors)[1]
        R = Diagonal(vec(fac[j,:])) * R
    end

    return GhostImage(gd.Proj, R, gd.flip_factors)
end


function extrapolate_values(x::Vector{Float64}, xb::Vector{Float64}, gd::GhostData)

    return gd.R*transpose(gc.Proj)*x + (I - gd.R)*xb
end


function require_boundary_conditions(A::SparseMatrixCSC, gd::GhostData)

    P = gd.Proj
    Pt = transpose(P)

    Anew = P*A*gd.R*Pt
    Bterm = P*A*(I - gd.R)

    return Anew, Bterm
end
