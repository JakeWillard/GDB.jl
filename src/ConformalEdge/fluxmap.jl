

function rk4_step(r::Vector{Float64}, f::Function, ds::Float64)

    k1 = f(r)
    k2 = f(r + ds*k1/2)
    k3 = f(r + ds*k2/2)
    k4 = f(r + ds*k3)

    return r + ds*(k1 + 2*k2 + 2*k3 + k4) / 6
end


struct FluxMap

    Z1 :: Matrix{Complex{Float64}}
    Z2 :: Matrix{Complex{Float64}}
    dS1 :: Matrix{Float64}
    dS2 :: Matrix{Float64}
    dr :: Float64
    dt :: Float64
end


function FluxMap(b::Function, Nr, Nt, deltaPhi, ds)

    rs = LinRange(0, 1, Nr)
    ts = LinRange(-pi, pi, Nt+1)[1:Nt]

    zs = pmap(1:Nr) do i
        col = zeros(Complex{Float64}, (4, Nt)
        for j=1:Nt
            x1 = rs[i]*[cos(ts[j]), sin(ts[j]), 0]
            x2 = x1[:]
            ds1 = 0.0
            ds2 = 0.0

            # trace forwards
            while x1[3] < deltaPhi
                x1 = rk4_step(x1, b, ds)
                ds1 += ds
            end

            # trace backwards
            while x2[3] < deltaPhi
                x2 = rk4_step(x2, b, -ds)
                ds2 += ds
            end

            # write final positions as complex numbers
            z1 = x1[1] + im*x1[2]
            z2 = x2[1] + im*x2[2]

            col[:,j] = [z1, z2, ds1, ds2]
        end
        col
    end
    Z1 = vcat([c[1:1,:] for c in zs]...)
    Z2 = vcat([c[2:2,:] for c in zs]...)
    dS1 = vcat([c[3:3,:] for c in zs]...)
    dS2 = vcat([c[4:4,:] for c in zs]...)

    return FluxMap(Z1, Z2, dS1, dS2, rs[2]-rs[1], ts[2]-ts[1])
end


function (fm::FluxMap)(z::Complex{Float64})

    r, t = abs(r), angle(z)
    i = Int64(floor(r / fm.dr))
    j = Int64(floor((t + pi) / fm.dt))
    u = r / fm.dr - i
    v = (t + pi) / fm.dt - j

    # bilinear interpolation
    z1 = [1 - u u] * fm.Z1[i:i+1,j:j+1] * [1 - v, v]
    z2 = [1 - u u] * fm.Z2[i:i+1,j:j+1] * [1 - v, v]
    ds1 = [1 - u u] * fm.dS1[i:i+1,j:j+1] * [1 - v, v]
    ds2 = [1 - u u] * fm.dS2[i:i+1,j:j+1] * [1 - v, v]
    return z1, z2, ds1, ds2
end
