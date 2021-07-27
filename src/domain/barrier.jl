

struct Barrier

    # surface defined as points where func(x,y) = val
    # rmap reflects points across the surface
    func::Function
    val::Float64
    rmap::Function

end


function smoothstep(x, y, delta, bar::Barrier)

    u = (bar.func(x,y) - bar.val) / delta + 0.5

    if u < 0
        return 0
    elseif 0 < u < 1
        return 3*u^2 - 2*u^3
    else
        return 1
    end
end


function reflection(x, y, delta, bar::Barrier)

    if 0 < smoothstep(x, y, delta, bar) < 1
        return bar.rmap(x, y)
    else
        return x, y
    end
end


function check_if_inside(x, y, deltas::Vector{Float64}, bars::Vector{Barrier})

    xy_is_inside = true

    for l=1:length(bars)
        if smoothstep(x, y, deltas[l], bars[l]) == 0
            xy_is_inside = false
            break
        end
    end

    return xy_is_inside
end


#convenience Barriers that will include all points or exclude all points
function Everywhere()

    func(x, y) = 0
    val = -Inf
    rmap(x, y) = x, y
    return Barrier(func, val, rmap)
end


function Nowhere()

    func(x, y) = 0
    val = Inf
    rmap(x, y) = x, y
    return Barrier(func, val, rmap)
end


# functions for geometric operations with polygons
function xy_to_uv(x, y, r1, r2)

    xs = x - r1[1]
    ys = y - r1[2]
    t = atan(r2[2] - r1[2], r2[1] - r1[1])

    u = xs * cos(t) + ys * sin(t)
    v = -xs * sin(t) + ys * cos(t)

    return Float64[u, v]
end


function uv_to_xy(u, v, r1, r2)

    t = atan(r2[2] - r1[2], r2[1] - r1[1])
    x = u * cos(t) - v * sin(t)
    y = u * sin(t) + v * cos(t)

    return r1 + Float64[x, y]
end


function poly_dist(x, y, vert)

    N = size(vert, 2) - 1
    h = zeros(Float64, N)

    for i=1:N
        r1 = vert[:,i]
        r2 = vert[:,i+1]
        u, v = xy_to_uv(x, y, r1, r2)
        if u < 0
            h[i] = norm([x,y] - r1) * sign(v)
        elseif 0 <= u < norm(r2 - r1)
            h[i] = v
        else
            h[i] = norm([x,y] - r2) * sign(v)
        end
    end

    h2, k = findmin(h .^ 2)
    return h[k], k
end


function ray_reflection(x, y, r1, r2)

    u, v = xy_to_uv(x, y, r1, r2)
    return uv_to_xy(u, -v, r1, r2)
end


# Barrier from polygon of given vertices
function PolyBarrier(vert)

    func(x,y) = poly_dist(x, y, vert)[1]

    function rmap(x, y)

        h, ns = poly_dist(x, y, vert)
        r1 = vert[:,ns]
        r2 = vert[:,ns+1]
        return ray_reflection(x, y, r1, r2)
    end

    return Wall(func, 0.0, rmap)
end
