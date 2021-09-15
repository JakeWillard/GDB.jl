




function shaped_transform(N, delta, kappa)

    epsilon = 1.0
    ts = LinRange(0, 2*pi, N+1)[1:N]
    alpha = asin(delta)
    x = epsilon * cos.(ts + alpha*sin.(ts))
    y = epsilon * kappa * sin.(ts)
    p = Polygon(x, y)
    m = ExteriorMap(p)

    f(x,y) = begin
        w = m(x + y*im)
        real(w), imag(w)
    end

    return f
end


function visualize_map(f::Function, L, N)

    N = 100
    xs = LinRange(-L/2, L/2, N)
    ys = LinRange(-L/2, L/2, N)
    z_real = zeros(N, N)
    z_imag = zeros(N, N)

    for i=1:N
        for j=1:N
            w = f(xs[i] + ys[j]*im)
            z_real[j, i] = real(w)
            z_imag[j, i] = imag(w)
        end
    end

    p = contour(xs, ys, z_real, levels=50) # levels=xs[1:skip:end], colorbar=false)
    contour!(xs, ys, z_imag, levels=50) # levels=ys[1:skip:end], colorbar=false)
    return p
end


function map_polar(f::Function, r, N)

    ts = LinRange(0, 2*pi, N)
    zs = r*(cos.(ts) + im*sin.(ts))
    ws = f.(zs)

    xs = real.(ws)
    ys = imag.(ws)
    return xs, ys
end


function map_d_shape(kappa, delta, Nr, Nt)

    alpha = asin(delta)
    rs = LinRange(1.0, 2.0, Nr)
    xs = zeros(Nt, Nr)
    ys = zeros(Nt, Nr)
    for i=1:Nr
        xs[:,i], ys[:,i] = map_polar(rs[i], Nt) do z
            elpart = (1 + kappa)*z/2.0 + (1 - kappa)/(2*z)
            sumpart = 0
            for j=1:5
                n = 2*j
                _a = (-1)^(n-2)*alpha^(n-1) / (factorial(n-1)*2^n)
                _b = 1 + (alpha/(2*n))*(z + 1/z)
                _c = (z - 1/z)^n
                sumpart += _a * _b * _c
            end
            elpart + sumpart
        end
    end

    return xs, ys
end
