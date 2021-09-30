
function generate_maps(ep, kapp, delta, Nv)

    alpha = asin(delta)

    ts = LinRange(0, 2*pi, Nv+1)[1:Nv]
    verts = zeros(Nv, 2)
    verts[:,1] = ep*cos.(ts + alpha*sin.(ts))
    verts[:,2] = ep*kapp*sin.(ts)

    f = ConformalMap(verts, 0.0)
    finv = inv(f)

    return f, finv
end


struct ExternalCoil

    Bref :: Float64
    refdist :: Float64
    r0 :: Vector{Float64}

end


function (c::ExternalCoil)(x::Float64, y::Float64)

    r = norm([x,y] - c.r0)
    return c.refdist*c.Bref*log(r / c.refdist)
end


function (c::ExternalCoil)(z::Complex{Float64})

    return c(real(z), imag(z))
end


function fourier_interp(dat)

    C = fft(dat)

    f(t) = begin

        Exps = [exp(2*pi*im*n*t) for n=1:50]
        Exps = [Exps; reverse(conj.(Exps))]
        return (C[1] + dot(C[2:end], Exps.^t)) / 101
    end

    return f


    return t -> ((C[1] + dot(C[2:end], Exps.^t)) / 101)
end


function generate_edge_field(Bp::Float64, coils::Vector{ExternalCoil}, ep, kapp, delta, Nv)

    f, finv = generate_maps(ep, kapp, delta, Nv)

    circle_boundary = [0.95*exp(im*t) for t in LinRange(0, 2*pi,102)[1:101]]
    shaped_boundary = finv.(circle_boundary)
    psi_b = Float64[sum([c(z) for c in coils]) for z in shaped_boundary]

    C = fft(psi_b) / 101

    function psi(x, y)
        z0 = x + im*y
        z = f(z0)
        V = [(z/0.95)^n for n=0:50]
        V = [V; 1 ./ reverse(V[1:50])]
        _a = real(dot(C, V))

        return _a + 0.95*Bp*log(abs(z)/0.95)
    end

    return psi
end
