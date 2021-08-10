

function example_geometry_setup(path::String, Nx, Ny)

    # define the magnetic field
    psi, bx, by, bz = solovev_flux_function(0.1, 0.1, 0.3, 1.7, 10.0, downsep=[1, -0.561])

    # values for flux surfaces
    psi_in = psi(1 - 0.28, 0)
    psi_out = psi(1 - 0.34, 0)
    psi_priv = psi(1, -(0.561 + 0.01))

    # define barriers for flux surfaces
    inner_flux = Barrier() do
        func(x,y) = (y > -0.561) ? psi(x,y) : Inf
        rmap(x,y) = trace_reflection(x, y, psi, psi_in, 0.001)
        func, psi_in, rmap, 1
    end

    outer_flux = Barrier() do
        func(x,y) = (norm([x-1,y+0.6]) > 0.05) ? psi(x,y) : psi_out + psi_priv - psi(x,y)
        rmap(x,y) = trace_reflection(x, y, func, psi_out, 0.001)
        func, psi_out, rmap, -1
    end

    # define barrier for divertor target plates
    target = Barrier() do
        func(x,y) = y
        val = -(0.561 + 0.04)
        rmap(x,y) = x, 2*val - y
        func, val, rmap, 1
    end

    # compute step function widths
    dr = 0.01
    dp_in = abs(ForwardDiff.derivative(u -> psi(u,0), 0.7)) * dr
    dp_out = abs(ForwardDiff.derivative(u -> psi(u,0), 0.6)) * dr
    dp_priv = abs(ForwardDiff.derivative(u -> psi(1,u), -(0.561 + 0.05))) * dr

    # create Grid
    deltas = Float64[dp_in, 0.5*(dp_out + dp_priv), dr]
    barrs = [inner_flux, outer_flux, target]
    grd = Grid(Nx, Ny, 1; Nbuffer=0) do
        inside(x,y) = check_if_inside(x, y, 2*deltas, barrs)
        r0 = Float64[0.5, -1]
        r1 = Float64[1.5, 1]
        inside, r0, r1
    end

    # compute penalization values as step functions
    pen = zeros(grd.Nk, 3)
    for i=1:3
        pen[:,i] = f_to_grid(grd) do x,y
            smoothstep(x, y, deltas[i], barrs[i])
        end
    end

    # K1 and K2 are operators for penalized forward-euler timesteps
    K1 = Diagonal(pen[:,1] .* pen[:,2] .* pen[:,3])
    K2 = 100 * (I - K1)

    # H1, H2, and H3 are for defining boundary values.
    h1 = h2 = h3 = ones(grd.Nk)
    h1[pen[:,1] .!= 1.0] .= 0.0
    h2[pen[:,2] .!= 1.0] .= 0.0
    h3[pen[:,3] .!= 1.0] .= 0.0
    H1 = Diagonal(h1)
    H2 = Diagonal(h2)
    H3 = Diagonal(h3)

    # compute ghost-conditions
    GC_dchlt = GhostConditions(2, 2, stencil2d(2, 2), barrs, grd)
    GC_nmann = swap_sign(GC_dchlt, 1, 2, 3)
    GC_u = swap_sign(GC_nmann, 3)

    # define vector for function that is +1 on left plate and -1 on right plate
    trgt_sgn = f_to_grid((x,y) -> sign(1 - x), grd)

    # lcfs average
    lcfs_avg = flux_surface_average(psi, 0.0, 150, 4, 4, stencil2d(4, 4), grd)

    # compute derivatives (need to renormalize since solovev assumes R=1, and we need a=1 consistent with the normalization factors in publications.)
    Minv = stencil1d(5)
    MinvT = stencil2d(4, 4)
    Dx = derivative_matrix(1, 0, Minv, Minv, grd) * 0.3
    Dy = derivative_matrix(0, 1, Minv, Minv, grd) * 0.3
    Dxy = derivative_matrix(1, 1, Minv, Minv, grd) * (0.3)^2
    Dxx = derivative_matrix(2, 0, Minv, Minv, grd) * (0.3)^2
    Dyy = derivative_matrix(0, 2, Minv, Minv, grd) * (0.3)^2
    Dxxx = derivative_matrix(3, 0, Minv, Minv, grd) * (0.3)^3
    Dyyy = derivative_matrix(0, 3, Minv, Minv, grd) * (0.3)^3
    Dxxy = derivative_matrix(2, 1, Minv, Minv, grd) * (0.3)^3
    Dxyy = derivative_matrix(1, 2, Minv, Minv, grd) * (0.3)^3
    Ds, Dss = fieldline_derivatives(bx, by, bz, 0.01, 4, 4, MinvT, 1, grd)

    return grd, GC_dchlt
end
