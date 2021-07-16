
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


struct Wall

    # surface defined as points where func(x,y) = val
    func::Function
    val::Float64
    reflect::Function

end


function PolyWall(vert)

    func(x,y) = poly_dist(x, y, vert)[1]

    function reflect(x, y)

        h, ns = poly_dist(x, y, vert)
        r1 = vert[:,ns]
        r2 = vert[:,ns+1]
        return ray_reflection(x, y, r1, r2)
    end

    return Wall(func, 0.0, reflect)
end


function FluxWall(psi, psi0, bx, by, ds)

    reflect(x, y) = trace_reflection(x, y, psi, bx, by, psi0, ds)

    return Wall(psi, psi0, reflect)
end


function smoothstep(x, y, delta, wall::Wall)

    u = (wall.func(x,y) - wall.val) / delta + 0.5

    if u < 0
        return 0
    elseif 0 < u < 1
        return 3*u^2 - 2*u^3
    else
        return 1
    end
end


function inside(x, y, deltas::Vector{Float64}, walls::Vector{Wall})

    is_inside = true

    for l=1:length(walls)
        if smoothstep(x, y, deltas[l], walls[l]) == 0
            is_inside = false
            break
        end
    end

    return is_inside
end
