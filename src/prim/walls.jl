
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


function smoothstep(x)

    if x < 0
        return 0
    elseif 0 < x < 1
        return 3*x^2 - 2*x^3
    else
        return 1
    end
end


struct Wall

    sstep::Function
    reflect::Function

end


function PolyWall(vert, ds)

    function sstep(x, y)

        h, ns = poly_dist(x, y, vert)
        u = h / ds
        return smoothstep(u)
    end

    function reflect(x, y)

        h, ns = poly_dist(x, y, vert)
        r1 = vert[:,ns]
        r2 = vert[:,ns+1]
        return ray_reflection(x, y, r1, r2)
    end

    return Wall(sstep, reflect)
end


function FluxWall(psi, psi0, bx, by, ds, dp)

    function sstep(x, y)

        u = (psi(x, y) - psi0) / dp
        return smoothstep(u)
    end

    reflect(x, y) = trace_reflection(x, y, psi, bx, by, psi0, ds)

    return Wall(sstep, reflect)
end
