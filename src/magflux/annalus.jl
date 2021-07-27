

function generate_annalus_field(a, R0, q0, s0, beta0)

    # define psi, q, and Bphi in terms of minor radius r
    psi(r) = ((r - a)/a + (1 - s0)/2 * ((r - a)/a)^2) / beta0
    q(r) = q0 / (1 - s0*(r - a)/a)

    # define in cylindrical coordiantes
    psi(x,y) = psi(sqrt((x-R0)^2 + y^2))
    q(x,y) = q(sqrt((x-R0)^2 + y^2))

    Bx(x, y) = ForwardDiff(u -> psi(x,u), y)
    By(x, y) = ForwardDiff(u -> -psi(u,y), x)
    Bphi(x, y) = (x / sqrt((x-R0)^2 + y^2)) * q(x,y) * norm(Float64[Bx(x,y)^2 + By(x,y)])
    B(x,y) = norm(Float64[Bx(x,y), By(x,y), Bphi(x,y)])

    # define unit vectors for b
    bx(x,y) = Bx(x,y) / B(x,y)
    by(x,y) = By(x,y) / B(x,y)
    bz(x,y) = Bz(x,y) / B(x,y)

    return psi, bx, by, bz
end


function AnnalusBarrier(x0, a, R0, q0, s0, beta0, ds)

    psi, bx, by, bz = generate_annalus_field(a, R0, q0, s0, beta0)
    psi0 = psi(x0, 0)
    rmap(x,y) = trace_reflection(x, y, psi, bx, by, psi0, ds)
    return Barrier(psi, psi0, rmap)
end


function LimiterBarrier(L, h, a, R0)

    # define vertices in counter-clockwise order
    vert = zeros(Float64, (2, 5))
    vert[:,1] = [R0-a-L, -h]
    vert[:,2] = [R0-a, -h]
    vert[:,3] = [R0-a, h]
    vert[:,4] = [R0-a-L, h]
    vert[:,4] = [R0-a-L, -h]
    return PolyBarrier(vert)
end
