
function rk4_step(r::Vector{Float64}, f::Function, ds::Float64)

    k1 = f(r)
    k2 = f(r + ds*k1/2)
    k3 = f(r + ds*k2/2)
    k4 = f(r + ds*k3)

    return r + ds*(k1 + 2*k2 + 2*k3 + k4) / 6
end

function trace(x0, y0, b::Function, deltaPhi, ds)

    x = [x0, y0, 0]
    deltaS = 0.0

    while x[3] < deltaPhi
        x[:] = rk4_step(x, b, ds)
        deltaS += abs(ds)
    end

    return x[1:2], deltaS
end


struct FluxMap

    Z :: Matrix{Complex{Float64}}
    dr :: Float64
    dt :: Float64
end


function FluxMap(b::Function, Nr, Nt, deltaPhi, ds)

    rs = LinRange(0, 1, Nr)
    ts = LinRange(-pi, pi, Nt+1)[1:Nt]
    Z = zeros(Nr, Nt)

    zs = pmap(1:Nr) do i
        col = zeros(Complex{Float64}, Nt)
        for j=1:Nt
            x, y = rs[i]*[cos(ts[j]), sin(ts[j])]
            xn, yn = trace(x, y, b, deltaPhi, ds)
            z[j] = xn + im*yn
        end
        col
    end
    Z = hcat(zs)

    return FluxMap(Z, rs[2]-rs[1], ts[2]-ts[1])
end


function (fm::FluxMap)(z::Complex{Float64})

    r, t = abs(r), angle(z)
    i = Int64(floor(r / fm.dr))
    j = Int64(floor(t / fm.dt))
    u = r / fm.dr - i
    v = t / fm.dt - j

    # bilinear interpolation
    zn = [1 - u u] * fm.Z[i:i+1,j:j+1] * [1 - v, v]
    return zn
end
