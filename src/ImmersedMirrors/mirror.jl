
function rk4_step(r::Vector{Float64}, f::Function, ds::Float64)

    k1 = f(r)
    k2 = f(r + ds*k1/2)
    k3 = f(r + ds*k2/2)
    k4 = f(r + ds*k3)

    return r + ds*(k1 + 2*k2 + 2*k3 + k4) / 6
end


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

    # find nearest arm, return distance to nearest arm
    dists = [distance_to_arm(x, y, M.arms[:,1,k], M.arms[:,2,k], M.arms[:,3,k]) for k=1:size(M.arms)[3]]
    _, kmin = findmin(dists .^ 2)

    return dists[kmin], kmin

    # # get the closest vertex, return distance to arm represented by that vertex
    # _, k = findmin([norm([x,y] - M.verts[:,i]) for i=1:size(M.verts)[2]])
    # vert0 = M.arms[:,1,k]
    # vert1 = M.arms[:,2,k]
    # vert2 = M.arms[:,3,k]
    #
    # return distance_to_arm(x, y, vert0, vert1, vert2), k
end


function mirror_image(x0, y0, dx, dy, M::Mirror)

    H(r) = distance_to_mirror(r[1], r[2], M)[1]
    gradH(r) = ForwardDiff.gradient(H, r)
    f(r) = gradH(r) / norm(gradH(r))
    x = x0
    y = y0
    ds = 0.1*minimum([dx, dy])
    dist = 0
    steps = 0

    while H([x,y]) < 0
        x, y = rk4_step([x,y], f, ds)
        dist += 1
        steps += 2
    end

    for _=1:dist
        x, y = rk4_step([x,y], f, ds)
    end

    i = round((x - x0)/dx)
    j = round((y - y0)/dy)
    while H([x0 + dx*i, y0 + dy*j]) < 0
        x, y = rk4_step([x,y], f, ds)
        i = round((x - x0)/dx)
        j = round((y - y0)/dy)
    end

    i = round((x - x0)/dx)
    j = round((y - y0)/dy)
    x = x0 + dx*i
    y = y0 + dy*j
    return x, y, i, j
end


function smoothstep(x, y, s0, ds, M::Mirror)

    u = (distance_to_mirror(x, y, M)[1] - s0) / ds
    return (u < 0) ? 0.0 : ((0 < u < 1) ? 3*u^2 - 2*u^3 : 1.0)
end


@recipe function f(M::Mirror)

    legend --> false
    seriescolor --> :black
    linewidth --> 3
    aspect_ratio --> :equal

    x = M.arms[1,:,:]
    y = M.arms[2,:,:]
    x, y
end
