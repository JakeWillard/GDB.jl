


function save_sparse_matrix(fid, A::SparseMatrixCSC, name::String)

    n, m = size(A)
    Is, Js, Vs = findnz(A)

    fid["$(name)_Is"] = Int32[Is...]
    fid["$(name)_Js"] = Int32[Js...]
    fid["$(name)_Vs"] = Float64[Vs...]
    fid["$(name)_size"] = Int32[n, m]

    @info "Sparse Matrix '$(name)' saved to disk."
end


function load_sparse_matrix(fid, name::String)

    if isempty(fid["$(name)_Is"])
        Is = Int32[1]
        Js = Int32[1]
        Vs = Float64[0]
    else
        Is = fid["$(name)_Is"][:]
        Js = fid["$(name)_Js"][:]
        Vs = fid["$(name)_Vs"][:]
    end
    n, m = fid["$(name)_size"][:]

    @info "Loaded Sparse Matrix '$(name)'."

    return sparse(Is, Js, Vs, n, m)
end


function save_grid(fid, grd::Grid, name::String)

    # create group for grid
    create_group(fid, name)

    # save projection
    save_sparse_matrix(fid, grd.Proj, "$(name)/Proj")
    fid["$(name)/r0"] = grd.r0[:]
    fid["$(name)/r1"] = grd.r1[:]
    fid["$(name)/points"] = grd.points[:,:]
    fid["$(name)/spacing"] = Float64[grd.dx, grd.dy]
    fid["$(name)/_size"] = Int64[grd.Nk, grd.Nz, grd._Nx, grd._Ny, grd._Nbuffer]
    fid["$(name)/_NaNs"] = grd._nan_outside_boundaries[:,:]

    @info "Grid '$(name)' saved to disk."

end


function load_grid(fid, name::String)

    Proj = load_sparse_matrix(fid, "$(name)/Proj")
    r0 = fid["$(name)/r0"][:]
    r1 = fid["$(name)/r1"][:]
    points = fid["$(name)/points"][:,:]
    dx, dy = fid["$(name)/spacing"][:]
    Nk, Nz, Nx, Ny, Nbuffer = fid["$(name)/_size"][:]
    NaNs = fid["$(name)/_NaNs"][:,:]

    @info "Loaded Grid '$(name)'."

    return Grid(r0, r1, points, Proj, dx, dy, Nk, Nz, Nx, Ny, Nbuffer, NaNs)
end


function save_ghost_conditions(fid, gc::GhostConditions, name::String)

    # create group
    create_group(fid, name)

    save_sparse_matrix(fid, gc.Proj, "$(name)/Proj")
    save_sparse_matrix(fid, gc.Mirror, "$(name)/Mirror")
    fid["$(name)/swaps"] = gc.swaps[:,:]
end


function load_ghost_conditions(fid, name::String)

    Proj = load_sparse_matrix(fid, "$(name)/Proj")
    Mirror = load_sparse_matrix(fid, "$(name)/Mirror")
    swaps = fid["$(name)/swaps"][:,:]

    return GhostConditions(Proj, Mirror, swaps)
end



function example_geometry_setup(path::String, Nx, Ny)

    # define the magnetic field
    psi, bx, by, bz = solovev_flux_function(0.1, 0.1, 0.3, 1.7, 10.0, downsep=[1, -0.561])

    # values for flux surfaces
    psi_in = psi(1 - 0.28, 0)
    psi_out = psi(1 - 0.34, 0)
    psi_priv = psi(1, -(0.561 + 0.01))

    # define barriers for flux surfaces
    inner_flux = Barrier() do
        func(x,y) = (y > -0.561) ? psi(x,y) : abs(psi(x,y))
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
    grd = Grid(Nx, Ny, 1) do
        inside(x,y) = check_if_inside(x, y, 10*deltas, barrs)
        r0 = Float64[0.5, -1]
        r1 = Float64[1.5, 1]
        inside, r0, r1
    end

    @info grd.Nk

    # compute penalization values as step functions
    pen = zeros(grd.Nk, 3)
    for i=1:3
        pen[:,i] = f_to_grid(grd) do x,y
            smoothstep(x, y, deltas[i], barrs[i])
        end
    end

    # K1 and K2 are operators for penalized forward-euler timesteps
    K1 = sparse(Diagonal(pen[:,1] .* pen[:,2] .* pen[:,3]))
    K2 = 100 * (I - K1)

    # H1, H2, and H3 are for defining boundary values.
    h1 = ones(grd.Nk)
    h2 = ones(grd.Nk)
    h3 = ones(grd.Nk)
    h1[pen[:,1] .!= 0.0] .= 0.0
    h2[pen[:,2] .!= 0.0] .= 0.0
    h3[pen[:,3] .!= 0.0] .= 0.0
    H1 = sparse(Diagonal(h1))
    H2 = sparse(Diagonal(h2))
    H3 = sparse(Diagonal(h3))

    # compute ghost-conditions
    GC_nmann = GhostConditions(barrs, grd) #XXX this is backwards now I think.
    GC_dchlt = swap_sign(GC_nmann, 1, 2, 3)
    GC_u = swap_sign(GC_nmann, 3)

    # define vector for function that is +1 on left plate and -1 on right plate
    trgt_sgn = f_to_grid((x,y) -> sign(1 - x), grd)

    # lcfs average
    lcfs_avg = flux_surface_average(psi, 0.0, 150, Float64[1, 0], 4, 4, stencil2d(4, 4), grd)

    # trace field-line images
    img = fieldline_images(bx, by, bz, 0.01, grd)

    fid = h5open(path, "w")
    save_grid(fid, grd, "Grid")
    save_sparse_matrix(fid, K1, "K1")
    save_sparse_matrix(fid, K2, "K2")
    save_sparse_matrix(fid, H1, "H1")
    save_sparse_matrix(fid, H2, "H2")
    save_sparse_matrix(fid, H3, "H3")
    save_ghost_conditions(fid, GC_dchlt, "GC_dchlt")
    save_ghost_conditions(fid, GC_nmann, "GC_nmann")
    save_ghost_conditions(fid, GC_u, "GC_u")
    fid["trgt_sgn"] = trgt_sgn[:]
    save_sparse_matrix(fid, lcfs_avg, "lcfs_avg")
    fid["FL_img"] = img[:,:]
    close(fid)
end


function example_physics_setup(R0, n0, T0, B0, init_path, geo_path)

    am, ad, ki, ke, er, eg, ev, de2, eta = dimensionless_parameters(0.3*R0, R0, n0, T0, B0)
    psi, bx, by, bz = solovev_flux_function(0.1, 0.1, 0.3, 1.7, 10.0, downsep=[1, -0.561])

    geofid = h5open(geo_path, "r")
    grd = load_grid(geofid, "Grid")
    GC_dchlt = load_ghost_conditions(geofid, "GC_dchlt")
    lcfs_avg = load_sparse_matrix(geofid, "lcfs_avg")
    H1 = load_sparse_matrix(geofid, "H1")
    H2 = load_sparse_matrix(geofid, "H2")
    H3 = load_sparse_matrix(geofid, "H3")
    FL_img = geofid["FL_img"][:,:]
    close(geofid)

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
    Ds, Dss = fieldline_derivatives(FL_img, 4, 4, MinvT, grd)

    psi_in = psi(1 - 0.28, 0)
    dp = abs(psi_in)

    n = Te = Ti = f_to_grid(grd) do x,y
        psi_xy = psi(x,y)
        (y > -0.561) ? (3/4)*(tanh((psi_in/2 - psi_xy)/dp)+1) + 0.5 : 0.5
    end
    uvec = f_to_grid(grd) do x,y
        0.5*(tanh((-0.4 - y)/0.05) + 1)*tanh((1-x)/0.05)
    end
    u = hcat(uvec, uvec, uvec)
    Sp = f_to_grid(grd) do x,y
        psi_xy = psi(x,y)
        (y > -0.561) ? 1/2 * (tanh((psi_in/2 - psi_xy)/dp) + 1) : 0.0
    end

    Pe = n .* Te
    Pi = n .* Ti

    lnn = lnTe = lnTi = begin v = log.(n); hcat(v, v, v) end
    j = jn = psi = zeros(grd.Nk)
    u = w = A = zeros(grd.Nk, 3)

    lnn_x = Dx * lnn[:,2]
    lnn_y = Dy * lnn[:,2]
    Pi_xx = Dxx * Pi
    Pi_yy = Dyy * Pi
    phi_b = bval_phi(w[:,2], Te, lcfs_avg, H1, H2, H3)
    phi = solve_vorticity_eqn(phi_b, w[:,2], n, lnn_x, lnn_y, Pi_xx, Pi_yy, ad, Dx, Dy, Dxx, Dyy, GC_dchlt)

    fid = h5open(init_path, "w")
    save_sparse_matrix(fid, Dx, "Dx")
    save_sparse_matrix(fid, Dy, "Dy")
    save_sparse_matrix(fid, Dxy, "Dxy")
    save_sparse_matrix(fid, Dxx, "Dxx")
    save_sparse_matrix(fid, Dyy, "Dyy")
    save_sparse_matrix(fid, Dxxx, "Dxxx")
    save_sparse_matrix(fid, Dyyy, "Dyyy")
    save_sparse_matrix(fid, Dxxy, "Dxxy")
    save_sparse_matrix(fid, Dxyy, "Dxyy")
    save_sparse_matrix(fid, Ds, "Ds")
    save_sparse_matrix(fid, Dss, "Dss")
    fid["lnn"] = lnn[:,:]
    fid["lnTe"] = lnTe[:,:]
    fid["lnTi"] = lnTi[:,:]
    fid["u"] = u[:,:]
    fid["w"] = w[:,:]
    fid["A"] = A[:,:]
    fid["phi"] = phi[:]
    fid["psi"] = psi[:]
    fid["n"] = n[:]
    fid["Te"] = Te[:]
    fid["Ti"] = Ti[:]
    fid["Pe"] = Pe[:]
    fid["Pi"] = Pi[:]
    fid["j"] = j[:]
    fid["jn"] = jn[:]
    fid["Sp"] = Sp[:]
    fid["params"] = Float64[am, ad, ki, ke, er, eg, ev, de2, eta]
    close(fid)
end


function test_simulate(Nt, sn, ste, sti, kdiff, dt, output_path, init_path, geo_path)

    geofid = h5open(geo_path, "r")
    grd = load_grid(geofid, "Grid")
    K1 = load_sparse_matrix(geofid, "K1")
    K2 = load_sparse_matrix(geofid, "K2")
    H1 = load_sparse_matrix(geofid, "H1")
    H2 = load_sparse_matrix(geofid, "H2")
    H3 = load_sparse_matrix(geofid, "H3")
    GC_dchlt = load_ghost_conditions(geofid, "GC_dchlt")
    GC_nmann = load_ghost_conditions(geofid, "GC_nmann")
    GC_u = load_ghost_conditions(geofid, "GC_u")
    trgt_sgn = geofid["trgt_sgn"][:]
    lcfs_avg = load_sparse_matrix(geofid, "lcfs_avg")
    FL_img = geofid["FL_img"][:,:]
    close(geofid)

    fid = h5open(init_path, "r")
    Dx = load_sparse_matrix(fid, "Dx")
    Dy = load_sparse_matrix(fid, "Dy")
    Dxy = load_sparse_matrix(fid, "Dxy")
    Dxx = load_sparse_matrix(fid, "Dxx")
    Dyy = load_sparse_matrix(fid, "Dyy")
    Dxxx = load_sparse_matrix(fid, "Dxxx")
    Dyyy = load_sparse_matrix(fid, "Dyyy")
    Dxxy = load_sparse_matrix(fid, "Dxxy")
    Dxyy = load_sparse_matrix(fid, "Dxyy")
    Ds = load_sparse_matrix(fid, "Ds")
    Dss = load_sparse_matrix(fid, "Dss")
    lnn = fid["lnn"][:,:]
    lnTe = fid["lnTe"][:,:]
    lnTi = fid["lnTi"][:,:]
    u = fid["u"][:,:]
    w = fid["w"][:,:]
    A = fid["A"][:,:]
    phi = fid["phi"][:]
    psi = fid["psi"][:]
    n = fid["n"][:]
    Te = fid["Te"][:]
    Ti = fid["Ti"][:]
    Pe = fid["Pe"][:]
    Pi = fid["Pi"][:]
    j = fid["j"][:]
    jn = fid["jn"][:]
    Sp = fid["Sp"][:]
    am, ad, ki, ke, er, eg, ev, de2, eta = fid["params"][:]
    close(fid)

    Sn = sn * Sp
    STe = ste * Sp
    STi = sti * Sp
    kdiff_lnn, kdiff_lnTe, kdiff_lnTi, kdiff_u, kdiff_w, kdiff_A = kdiff[:]
    return phi, grd
    for t=1:Nt
        leapfrog!(1, lnn, lnTe, lnTi, u, w, A, phi, psi, n, Te, Ti, Pe, Pi, j, jn, Sn, STe, STi, Dx, Dy, Dxy, Dxx, Dyy, Dxxx, Dyyy,
                  Dxxy, Dxyy, Ds, Dss, dt, 10, K1, K2, H1, H2, H3, trgt_sgn, lcfs_avg, GC_nmann, GC_dchlt, GC_u, kdiff_lnn,
                  kdiff_lnTe, kdiff_lnTi, kdiff_u, kdiff_w, kdiff_A, am, ad, ki, ke, er, eg, ev, de2, eta)
        # write to output
        # ;alksdjfl;kasjfl;sk
    end

    return n, grd
end
