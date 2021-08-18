
function find_surface(r0, dr, psi, psi0)

    r = r0[:]
    sgn = psi(r0...) - psi0
    while true
        rnew += dr
        if sgn*(psi(rnew...) - psi0) < 0
            return 0.5*(rnew + r)
        end
        r[:] = rnew
    end
end


function regula_falsi(f::Function, val::Float64, bracket1::Vector{Float64}, bracket2::Vector{Float64})

    a = bracket1[:]
    b = bracket2[:]
    fa = f(a) - val
    fb = f(b) - val
    c = (a*fb - b*fa) / (fb - fa)
    fc = f(c) - val

    while abs(fc) > 1e-8
        if fa*fc < 0
            b = c
        elseif fc*fb < 0
            a = c
        else
            @error "root finding algorithm has failed."
            break
        end
        fa = f(a) - val
        fb = f(b) - val
        c = (a*fb - b*fa) / (fb - fa)
        fc = f(c) - val
    end

    return c
end


function rk4_step(r::Vector{Float64}, f::Function, ds::Float64)

    k1 = f(r)
    k2 = f(r + ds*k1/2)
    k3 = f(r + ds*k2/2)
    k4 = f(r + ds*k3)

    return r + ds*(k1 + 2*k2 + 2*k3 + k4) / 6
end


function conservative_step(r::Vector{Float64}, b::Function, psi::Function, ds::Float64)

    # compute predictor step using rk4
    rp = rk4_step(r, b, ds)

    # compute correction using root finding algorithm
    en = [0 -1; 1 0] * b(rp)
    bracket1 = rp + ds*en
    bracket2 = rp - ds*en
    rc = regula_falsi(psi, psi(r), bracket1, bracket2)

    # compute with correct displacement
    return rc, norm(rc - r)
end



function trace_fieldline(stop_condition::Function, r0::Vector{Float64}, b::Function, psi::Function, ds::Float64)

    r = r0[:]
    deltaS = 0.0

    points = Vector{Float64}[r0]
    while !stop_condition(r)
        rnew, dsnew = conservative_step(r, b, psi, ds)
        r[:] = rnew
        deltaS += dsnew
        append!(points, [rnew])
    end

    return hcat(points...), deltaS
end


function trace_fieldline(N::Int64, r0::Vector{Float64}, b::Function, psi::Function, ds::Float64)

    r = r0[:]
    deltaS = 0.0

    points = zeros(2, N)
    points[:,1] = r[:]
    for s=2:N
        rnew, dsnew = conservative_step(r, b, psi, ds)
        r[:] = rnew
        deltaS += dsnew
        points[:,s] = r
    end

    return points, deltaS
end
