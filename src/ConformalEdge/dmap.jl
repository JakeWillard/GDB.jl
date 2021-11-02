struct DMap

    params :: Tuple{Float64, Float64, Float64}
    d_to_c :: ConformalMap{Float64}
    c_to_d :: InverseConformalMap{Float64}

end


function DMap(ep, kapp, delta, Nv)

    params = (ep, kapp, delta)
    alpha = asin(delta)

    ts = LinRange(0, 2*pi, Nv+1)[1:Nv]
    verts = zeros(Nv, 2)
    verts[:,1] = ep*cos.(ts + alpha*sin.(ts))
    verts[:,2] = ep*kapp*sin.(ts)

    d_to_c = ConformalMap(verts, 0.0)
    c_to_d = inv(d_to_c)

    return DMap(params, d_to_c, c_to_d)
end


function (dm::DMap)(x::Float64, y::Float64)
    return dm.d_to_c(x + im*y)
end


function (dm::DMap)(z::Complex{Float64})
    zn = dm.c_to_d(z)
    return real(zn), imag(zn)
end


# function (dm::DMap)(psib::Function)
#
#     rs = dm.([exp(im*t) for t in LinRange(0, 2*pi, 1002)[1:1001]])
#     ts = [atan(r[2], r[1]) for r in rs]
#     psib_vec = psib.(ts)
#     C = fft(psib_vec)
#
#     c0 = real(C[1]) / 1001
#     a = 2*real.(C[2:501]) / 1001
#     b = -2*imag.(C[2:501]) / 1001
#
#     psi_v(z) = begin
#         r = abs(z)
#         t = angle(z)
#         return c0 + sum([(r/0.95)^k * (a[k]*cos(k*t) + b[k]*sin(k*t)) for k=1:500])
#     end
# end
