

function annalus_long_calculations(L, h, a, R0, q0, s0, beta0, ds, delta, N, Nz, m, path)

    # compute grid

    inner_flux_surface = AnnalusBarrier(R0-a+L/2, a, R0, q0, s0, beta0, ds)
    outer_flux_surface = AnnalusBarrier(R0-a-L/2, a, R0, q0, s0, beta0, ds)
    limiter_surface = LimiterBarrier(L, h, a, R0)

    deltas = Float64[delta, delta, delta]
    bars = [inner_flux_surface, outer_flux_surface, limiter_surface]
    inside(x, y) = check_if_inside(x, y, deltas, bars)
    r0 = Float64[R0 - a - L/2 - 3*delta, R0 - a - L/2 - 3*delta]
    r1 = Float64[R0 + a + L/2 + 3*delta, R0 + a + L/2 + 3*delta]

    grd = Grid(inside, r0, r1, N, N)

    # compute fieldline derivatives

    psi, bx, by, bz = generate_annalus_fields(a, R0, q0, s0, beta0)
    MinvT = stencil2d(m ,m)
    Ds, Dss = fieldline_derivatives(bx, by, bz, ds, m, m, MinvT, Nz, grd)

    fid = h5open(path, "r+")
    save_grid(fid, "Grid", grd)
    save_sparse_matrix(fid, Ds, "Ds")
    save_sparse_matrix(fid, Dss, "Dss")
    close(fid)

end
