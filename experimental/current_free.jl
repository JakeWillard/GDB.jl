

function FluxBarrier(psi, psi0, bx, by, ds)
    rmap(x,y) = trace_reflection(x, y, psi, bx, by, psi0, ds)
    return Barrier(psi, psi0, rmap)
end


function current_free_barriers(Nv, L, epsilon, delta, kappa, ds)

    # ds: small distance for tracing (step size for rk4 integration)
    # Nv: number of vertices for outer wall
    # L: thickness of domain on inboard side

    alpha = asin(delta)
    x(theta) = 1 + epsilon*cos(theta + alpha*sin(theta))
    y(theta) = epsilon*kappa*sin(theta)

    # generate solovev field
    psi, bx, by, bz = generate_solovev_field()
    psi_in = psi(1 - epsilon, 0)
    psi_out = psi(1 - epsilon - L, 0)

    # inner and outer barriers
    inner_flux_surface = FluxBarrier(psi, psi_in, bx, by, ds)
    outer_flux_surface = FluxBarrier(psi, psi_out, bx, by, ds)

    # outer wall
    vert = zeros(Float64, (2, Nv))
    thetas = LinRange(0, -2*pi, Nv)
    for i=1:Nv
        vert[1,:] = x(thetas[i])
        vert[2,:] = y(thetas[i])
    end
    outer_wall = PolyBarrier(vert)

    return Barrier[inner_flux_surace, outer_flux_surface, outer_flux_surface]
end
