

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
    r0 :: Float64
    dr :: Float64
    dt :: Float64
    deltaPhi :: Float64
end


function FluxMap(b::Function, r0, Nr, Nt, deltaPhi, ds)

    rs = LinRange(r0, 1, Nr)
    ts = LinRange(-pi, pi, Nt+1)[1:Nt]

    chnl = RemoteChannel(()->Channel{Bool}(), 1)
    p = Progress(Nr*Nt, desc="Tracing fieldlines... ")
    @async while take!(chnl)
        next!(p)
    end

    zs = pmap(1:Nr) do i
        col = zeros(Complex{Float64}, (4, Nt))
        for j=1:Nt
            x1 = rs[i]*[cos(ts[j]), sin(ts[j]), 0]
            x2 = x1[:]
            ds1 = 0.0
            ds2 = 0.0

            # trace forwards
            while abs(x1[3]) < deltaPhi
                x1 = rk4_step(x1, b, ds)
                ds1 += ds
            end

            # trace backwards
            while abs(x2[3]) < deltaPhi
                x2 = rk4_step(x2, b, -ds)
                ds2 += ds
            end

            # write final positions as complex numbers
            z1 = x1[1] + im*x1[2]
            z2 = x2[1] + im*x2[2]

            put!(chnl, true)
            col[:,j] = [z1, z2, ds1, ds2]
        end
        col
    end
    Z1 = vcat([c[1:1,:] for c in zs]...)
    Z2 = vcat([c[2:2,:] for c in zs]...)
    dS1 = vcat([c[3:3,:] for c in zs]...)
    dS2 = vcat([c[4:4,:] for c in zs]...)

    # make periodic
    Z1 = hcat(Z1, Z1[:,1])
    Z2 = hcat(Z2, Z2[:,1])
    dS1 = hcat(dS1, dS1[:,1])
    dS2 = hcat(dS2, dS2[:,1])

    return FluxMap(Z1, Z2, dS1, dS2, r0, rs[2]-rs[1], ts[2]-ts[1], deltaPhi)
end


function (fm::FluxMap)(z::Complex{Float64})

    r, t = abs(z), angle(z)
    i = Int64(floor((r - fm.r0) / fm.dr)) + 1
    j = Int64(floor(((t + pi) / fm.dt))) + 1
    u = (r - fm.r0) / fm.dr - i + 1
    v = (t + pi) / fm.dt - j + 1

    if (i > size(fm.Z1)[1]-1) || (j > size(fm.Z2)[2]-1) || (i < 1)
        return z, z, fm.deltaPhi, fm.deltaPhi
    end

    # bilinear interpolation
    jmax = (j==size(fm.Z1)[2]) ? 1 : j + 1
    z1 = dot([1-u, u], fm.Z1[i:i+1,[j, jmax]] * [1 - v, v])
    z2 = dot([1-u, u], fm.Z2[i:i+1,[j, jmax]] * [1 - v, v])
    ds1 = dot([1-u, u], fm.dS1[i:i+1,[j, jmax]] * [1 - v, v])
    ds2 = dot([1-u, u], fm.dS2[i:i+1,[j, jmax]] * [1 - v, v])
    return z1, z2, ds1, ds2
end


# function plot_arrows(fm::FluxMap, Nr::Int64, Nt::Int64)
#
#     rs = LinRange(fm.r0, 0.8, Nr)
#     ts = LinRange(-pi, pi, Nt+1)[1:Nt]
#
#     p = plot()
#     for r in rs
#         for t in ts
#             z1 = r*exp(im*t)
#             z2, z0 = fm(z1)[1:2]
#             plot!(p, [z0, z1, z2], arrow=true, legend=false, color=:black, linewidth=0.2)
#         end
#     end
#
#     return p
# end
