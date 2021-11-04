
# struct Hat
#
#     x0 :: Float64
#     x1 :: Float64
#     v :: Float64
# end
#
#
# function (h::Hat)(x::Float64)
#     u = abs(2*(x - (h.x0 + h.x1)/2) / (h.x1 - h.x0))
#     u < 1 ? h.v*(1 - u) : 0
# end


struct SOL

    inner_boundary :: Vector{Complex{Float64}}
    outer_boundary :: Vector{Complex{Float64}}
    psi :: Function
    qinv :: Function
    fm :: FluxMap

end


function SOL(Bpol, Bz, Nc, r0, Nr, Nt, Nz, ds, divs::Matrix{Float64})

    # compute inner and outer boundaries
    outer_chain = Complex{Float64}[]
    inner_chain = Vector{Complex{Float64}}(undef, Nc)
    thetas = LinRange(-pi, pi, Nc+1)[1:Nc]
    for i=1:Nc
        if !any([divs[1,j] < thetas[i] < divs[2,j] for j=1:size(divs)[2]])
            outer_chain = [outer_chain; [0.95*exp(im*thetas[i])]]
        end
        inner_chain[i] = r0*exp(im*thetas[i])
    end

    # compute boundary value for psi
    psib = Vector{Float64}(undef, 1001)
    ts = LinRange(0, 2*pi, 1002)[1:1001]
    for i=1:1001
        for j=1:size(divs)[2]
            if divs[1,j] < ts[i] < divs[2,j]
                u = abs(2*(ts[i] - (divs[1,j] + divs[2,j])/2) / (divs[2,j] - divs[1,j]))
                psib[i] = divs[3,j]*(1 - u)
            else
                psib[i] = 0.0
            end
        end
    end

    # find fourier coefficients for flux function
    C = fft(psib)
    c0 = real(C[1]) / 1001
    a = 2*real.(C[2:501]) / 1001
    b = -2*imag.(C[2:501]) / 1001

    # define flux function
    psi(z) = begin
        r = abs(z)
        t = angle(z)

        kmax = minimum([500, abs(Int64(ceil(8 / log10(r/0.9))))])
        if r < 0.9
            c0 + Bpol*log(r/0.9) + sum([(r/0.9)^k * (a[k]*cos(k*t) + b[k]*sin(k*t)) for k=1:kmax])
        else
            c0 + Bpol*log(r/0.9) + sum([(r/0.9)^(-k) * (a[k]*cos(k*t) + b[k]*sin(k*t)) for k=1:kmax])
        end
    end

    # define poloidal and total magnetic field
    gradPsi(r) = ForwardDiff.gradient(r -> psi(r[1] + im*r[2]), r)
    Bp(r) = begin dpsi = gradPsi(r); [-dpsi[2], dpsi[1]] end
    B(r) = [Bp(r)..., Bz]
    bvec(r) = begin Bvec = B(r); Bvec / norm(Bvec) end

    # define 1/q
    qinv(z) = norm(Bp([real(z), imag(z)])) / Bz

    # compute flux map
    fm = FluxMap(bvec, 0.95*r0, Nr, Nt, 2*pi/Nz, ds)

    return SOL(inner_chain, outer_chain, psi, qinv, fm)
end
