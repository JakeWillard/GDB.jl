

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


function save_ghost_data(fid, gd::GhostData, name::String)

    # create group
    create_group(fid, name)

    save_sparse_matrix(fid, gd.Proj, "$(name)/Proj")
    save_sparse_matrix(fid, gd.Mirror, "$(name)/Mirror")
    fid["$(name)/flip"] = gd.flip_factors[:,:]
end


function load_ghost_data(fid, name::String)

    Proj = load_sparse_matrix(fid, "$(name)/Proj")
    Mirror = load_sparse_matrix(fid, "$(name)/Mirror")
    flip_factors = fid["$(name)/flip"][:,:]

    return GhostData(Proj, Mirror, flip_factors)
end


function boundary_setup(bdry_path::String, N1, N2, N3)

    count = 0

    # define the magnetic field
    psi, bx, by, bz = solovev_flux_function(0.1, 0.1, 0.3, 1.7, 10.0, downsep=[1, -0.561])
    bp(r) = Float64[bx(r...), by(r...)]

    # values for flux surfaces
    psi_in = psi(1 - 0.28, 0)
    psi_out = psi(1 - 0.34, 0)
    psi_priv = psi(1, -(0.561 + 0.01))

    # find where the inner flux surface intersects y=0 on the outboard side
    bracket1 = [1.0, 0.0]
    bracket2 = [2.0, 0.0]
    r0 = rootfind_bisection(r->psi(r...), psi_in, bracket1, bracket2)

    # resolve inner boundary
    inner_chain_1, _ = trace_fieldline(r0, bp, r->psi(r...), 0.001) do r
        r[2] < 0
    end
    inner_chain_2, _ = trace_fieldline(inner_chain_1[:,end], bp, r->psi(r...), 0.001) do r
        r[2] > 0
    end
    skip = maximum([div(size(inner_chain_1)[2], N1), 1])
    inner_chain = hcat(inner_chain_1[:,1:skip:end-1], inner_chain_2[:,1:skip:end-1], inner_chain_1[:,1])
    Nin = size(inner_chain)[2] - 1

    inner_range = 1:Nin
    count += Nin

    # find where outer flux surface intersects right diverter plate
    bracket1 = [1.0, -(0.561 + 0.04)]
    bracket2 = [2.0, -(0.561 + 0.04)]
    r0 = rootfind_bisection(r->psi(r...), psi_out, bracket1, bracket2)

    # resolve outer flux surface
    outer_chain, _ = trace_fieldline(r0, bp, r->psi(r...), 0.001) do r
        r[2] < -(0.561 + 0.04)
    end
    skip = maximum([div(size(outer_chain)[2], N2), 1])
    outer_chain = outer_chain[:,1:skip:end]
    Nout = size(outer_chain)[2]

    outer_range_1 = count+1:count+Nout
    count += Nout

    # find where outer flux surface intersects left diverter plate, add to outer_chain
    bracket1 = [0.0, -(0.561 + 0.04)]
    bracket2 = [1.0, -(0.561 + 0.04)]
    rl = rootfind_bisection(r->psi(r...), psi_out, bracket1, bracket2)
    outer_chain = hcat(outer_chain, rl)

    div_range_1 = count+1:count+1
    count += 1

    # find where private flux surface intersects left diverter plate
    bracket1 = [rl[1], -(0.561 + 0.04)]
    bracket2 = [1.0, -(0.561 + 0.04)]
    rl = rootfind_bisection(r->psi(r...), psi_priv, bracket1, bracket2)

    # trace private flux surface
    private_chain, _ = trace_fieldline(rl, r->-bp(r), r->psi(r...), 0.001) do r
        r[2] < -(0.561 + 0.04)
    end
    skip = maximum([div(size(private_chain)[2], N3), 1])
    private_chain = private_chain[:,1:skip:end]
    outer_chain = hcat(outer_chain, private_chain)
    Npriv = size(private_chain)[2]

    outer_range_2 = count+1:count+Npriv
    count += Npriv

    # find where private flux surface intersects right diverter plate
    bracket1 = [1.0, -(0.561 + 0.04)]
    bracket2 = [2.0, -(0.561 + 0.04)]
    rr = rootfind_bisection(r->psi(r...), psi_priv, bracket1, bracket2)
    outer_chain = hcat(outer_chain, rr, outer_chain[:,1])

    div_range_2 = count+1:count+1
    count += 1

    Nv = Nin + Nout + Npriv + 2

    inner_boundary = fill(false, Nv)
    inner_boundary[inner_range] .= true
    outer_boundary = fill(false, Nv)
    outer_boundary[outer_range_1] .= true
    outer_boundary[outer_range_2] .= true
    div_boundary = fill(false, Nv)
    div_boundary[div_range_1] .= true
    div_boundary[div_range_2] .= true

    verts = zeros(2, 2, Nv)
    verts[:,1,1:Nin] = inner_chain[:,1:Nin]
    verts[:,2,1:Nin] = inner_chain[:,2:Nin+1]
    verts[:,1,Nin+1:end] = outer_chain[:,1:end-1]
    verts[:,2,Nin+1:end] = outer_chain[:,2:end]
    orientations = ones(Int64, Nv)
    orientations[1:Nin] .= -1

    fid = h5open(bdry_path, "w")
    fid["verts"] = verts[:,:,:]
    fid["orientations"] = orientations[:]
    fid["inner_boundary"] = inner_boundary[:]
    fid["outer_boundary"] = outer_boundary[:]
    fid["div_boundary"] = div_boundary[:]
    close(fid)
end


function simple_geometry_setup(geo_path::String, bdry_path::String, Nx, Ny)

    fid = h5open(bdry_path, "r")
    verts = fid["verts"][:,:,:]
    orientations = fid["orientations"][:]
    inner_boundary = fid["inner_boundary"][:]
    outer_boundary = fid["outer_boundary"][:]
    div_boundary = fid["div_boundary"][:]
    close(fid)

    ms = MirrorSegments(verts, orientations)

    grd = Grid(Nx, Ny, 1) do
        inside(x,y) = distance_to_mirror(x, y, ms)[1] > -0.05
        r0 = Float64[0.5, -1]
        r1 = Float64[1.5, 1]
        inside, r0, r1
    end

    # compute penalization values as step functions
    pen = zeros(grd.Nk, 3)
    pen[:,1] = f_to_grid(grd) do x,y
        smoothstep(x, y, 0.02, 0.02, ms[inner_boundary])
    end
    pen[:,2] = f_to_grid(grd) do x,y
        smoothstep(x, y, 0.02, 0.02, ms[outer_boundary])
    end
    pen[:,3] = f_to_grid(grd) do x,y
        smoothstep(x, y, 0.02, 0.02, ms[div_boundary])
    end

    # K1 and K2 are operators for penalized forward-euler timesteps
    K1 = sparse(Diagonal(pen[:,1] .* pen[:,2] .* pen[:,3]))
    K2 = 100 * (I - K1)

    gd_nmann = GhostData(ms, grd)

    fid = h5open(geo_path, "w")
    save_grid(fid, grd, "Grid")
    save_sparse_matrix(fid, K1, "K1")
    save_sparse_matrix(fid, K2, "K2")
    save_ghost_data(fid, gd_nmann, "gd")
    close(fid)
end


function simple_physics_setup(w0, psi0, L_w, L_psi, init_path::String, geo_path::String)

    fid = h5open(geo_path, "r")
    grd = load_grid(fid, "Grid")
    gd = load_ghost_data(fid, "gd")
    close(fid)

    psi_func, bx, by, bz = solovev_flux_function(0.1, 0.1, 0.3, 1.7, 10.0, downsep=[1, -0.561])
    psi_mid = psi_func(0.65, 0)
    dp_w = abs(ForwardDiff.derivative(u -> psi_func(u,0), 0.65)) * L_w
    dp_psi = dp_w * L_psi / L_w

    w_vec = f_to_grid(grd) do x,y
        w0*tanh((psi_func(x,y) - psi_mid) / dp_w)
    end
    psi_vec = f_to_grid(grd) do x,y
        psi0*tanh((psi_func(x,y) - psi_mid) / dp_psi)
    end

    w = hcat(w_vec, w_vec, w_vec)
    psi = hcat(psi_vec, psi_vec, psi_vec)

    # compute derivatives
    Minv = stencil1d(3)
    Dx = derivative_matrix(1, 0, Minv, Minv, grd) * 0.3
    Dy = derivative_matrix(0, 1, Minv, Minv, grd) * 0.3
    Dxx = derivative_matrix(2, 0, Minv, Minv, grd) * (0.3)^2
    Dyy = derivative_matrix(0, 2, Minv, Minv, grd) * (0.3)^2

    # solve vorticity equation
    phi = solve_vorticity_eqn(zeros(grd.Nk), w[:,2], Dxx, Dyy, gd)

    fid = h5open(init_path, "w")
    save_sparse_matrix(fid, Dx, "Dx")
    save_sparse_matrix(fid, Dy, "Dy")
    save_sparse_matrix(fid, Dxx, "Dxx")
    save_sparse_matrix(fid, Dyy, "Dyy")
    fid["w"] = w[:,:]
    fid["phi"] = phi[:]
    fid["psi"] = psi[:,:]
    close(fid)
end


function simulate(Nt, k_w, k_psi, am, dt, output_path, init_path, geo_path)

    fid = h5open(geo_path, "r")
    grd = load_grid(fid, "Grid")
    K1 = load_sparse_matrix(fid, "K1")
    K2 = load_sparse_matrix(fid, "K2")
    gd = load_ghost_data(fid, "gd")
    close(fid)

    fid = h5open(init_path, "r")
    Dx = load_sparse_matrix(fid, "Dx")
    Dy = load_sparse_matrix(fid, "Dy")
    Dxx = load_sparse_matrix(fid, "Dxx")
    Dyy = load_sparse_matrix(fid, "Dyy")
    w = fid["w"][:,:]
    phi = fid["phi"][:]
    psi = fid["psi"][:,:]
    close(fid)

    fid = h5open(output_path, "w")
    save_grid(fid, grd, "Grid")
    create_dataset(fid, "w", datatype(Float64), dataspace(grd.Nk*grd.Nz,Nt), chunk=(grd.Nk*grd.Nz,1))
    create_dataset(fid, "phi", datatype(Float64), dataspace(grd.Nk*grd.Nz,Nt), chunk=(grd.Nk*grd.Nz,1))
    create_dataset(fid, "psi", datatype(Float64), dataspace(grd.Nk*grd.Nz,Nt), chunk=(grd.Nk*grd.Nz,1))
    fid["w"][:,1] = w[:,2]
    fid["phi"][:,1] = phi[:]
    fid["psi"][:,1] = psi[:,2]
    close(fid)

    Nskip = 1
    p = Progress(Nskip*(Nt-1), "Running simple model... ")
    for t=2:Nt
        for _=1:Nskip
            leapfrog!(w, psi, phi, Dx, Dxx, Dy, Dyy, K1, K2, gd, dt, am, k_psi, k_w)
            next!(p)
        end

        fid = h5open(output_path, "r+")
        fid["w"][:,t] = w[:,2]
        fid["phi"][:,t] = phi[:]
        fid["psi"][:,t] = psi[:,2]
        close(fid)
    end

end






#
#
#
# function simple_geometry_setup(path::String, Nx, Ny)
#
#     # define the magnetic field
#     psi, bx, by, bz = solovev_flux_function(0.1, 0.1, 0.3, 1.7, 10.0, downsep=[1, -0.561])
#
#     # values for flux surfaces
#     psi_in = psi(1 - 0.28, 0)
#     psi_out = psi(1 - 0.34, 0)
#     psi_priv = psi(1, -(0.561 + 0.01))
#
#     # define barriers for flux surfaces
#     inner_flux = Barrier() do
#         func(x,y) = (y > -0.561) ? psi(x,y) : abs(psi(x,y))
#         rmap(x,y) = trace_reflection(x, y, psi, psi_in, 0.001)
#         func, psi_in, rmap, 1
#     end
#
#     outer_flux = Barrier() do
#         func(x,y) = (norm([x-1,y+0.6]) > 0.05) ? psi(x,y) : psi_out + psi_priv - psi(x,y)
#         rmap(x,y) = trace_reflection(x, y, func, psi_out, 0.001)
#         func, psi_out, rmap, -1
#     end
#
#     # define barrier for divertor target plates
#     target = Barrier() do
#         func(x,y) = y
#         val = -(0.561 + 0.04)
#         rmap(x,y) = x, 2*val - y
#         func, val, rmap, 1
#     end
#
#     # compute step function widths
#     dr = 0.01
#     dp_in = abs(ForwardDiff.derivative(u -> psi(u,0), 0.7)) * dr
#     dp_out = abs(ForwardDiff.derivative(u -> psi(u,0), 0.6)) * dr
#     dp_priv = abs(ForwardDiff.derivative(u -> psi(1,u), -(0.561 + 0.05))) * dr
#
#     # create Grid
#     deltas = Float64[dp_in, 0.5*(dp_out + dp_priv), dr]
#     barrs = [inner_flux, outer_flux, target]
#     grd = Grid(Nx, Ny, 1) do
#         inside(x,y) = check_if_inside(x, y, 10*deltas, barrs)
#         r0 = Float64[0.5, -1]
#         r1 = Float64[1.5, 1]
#         inside, r0, r1
#     end
#
#     # compute penalization values as step functions
#     pen = zeros(grd.Nk, 3)
#     for i=1:3
#         pen[:,i] = f_to_grid(grd) do x,y
#             smoothstep(x, y, deltas[i], barrs[i])
#         end
#     end
#
#     # K1 and K2 are operators for penalized forward-euler timesteps
#     K1 = sparse(Diagonal(pen[:,1] .* pen[:,2] .* pen[:,3]))
#     K2 = 100 * (I - K1)
#
#     gc_nmann = GhostConditions(barrs, grd)
#     gc = swap_sign(gc_nmann, 1, 2, 3)
#
#     fid = h5open(path, "w")
#     save_grid(fid, grd, "Grid")
#     save_sparse_matrix(fid, K1, "K1")
#     save_sparse_matrix(fid, K2, "K2")
#     save_ghost_conditions(fid, gc_nmann, "gc")
#     close(fid)
# end
#
#
# function simple_physics_setup(w0, psi0, L_w, L_psi, init_path::String, geo_path::String)
#
#     fid = h5open(geo_path, "r")
#     grd = load_grid(fid, "Grid")
#     gc = load_ghost_conditions(fid, "gc")
#     close(fid)
#
#     psi_func, bx, by, bz = solovev_flux_function(0.1, 0.1, 0.3, 1.7, 10.0, downsep=[1, -0.561])
#     psi_mid = psi_func(0.65, 0)
#     dp_w = abs(ForwardDiff.derivative(u -> psi_func(u,0), 0.65)) * L_w
#     dp_psi = dp_w * L_psi / L_w
#
#     w_vec = f_to_grid(grd) do x,y
#         w0*tanh((psi_func(x,y) - psi_mid) / dp_w)
#     end
#     psi_vec = f_to_grid(grd) do x,y
#         psi0*tanh((psi_func(x,y) - psi_mid) / dp_psi)
#     end
#
#     w = hcat(w_vec, w_vec, w_vec)
#     psi = hcat(psi_vec, psi_vec, psi_vec)
#
#     # compute derivatives
#     Minv = stencil1d(3)
#     Dx = derivative_matrix(1, 0, Minv, Minv, grd) * 0.3
#     Dy = derivative_matrix(0, 1, Minv, Minv, grd) * 0.3
#     Dxx = derivative_matrix(2, 0, Minv, Minv, grd) * (0.3)^2
#     Dyy = derivative_matrix(0, 2, Minv, Minv, grd) * (0.3)^2
#
#     # solve vorticity equation
#     phi = solve_vorticity_eqn(zeros(grd.Nk), w[:,2], Dxx, Dyy, gc)
#
#     fid = h5open(init_path, "w")
#     save_sparse_matrix(fid, Dx, "Dx")
#     save_sparse_matrix(fid, Dy, "Dy")
#     save_sparse_matrix(fid, Dxx, "Dxx")
#     save_sparse_matrix(fid, Dyy, "Dyy")
#     fid["w"] = w[:,:]
#     fid["phi"] = phi[:]
#     fid["psi"] = psi[:,:]
#     close(fid)
# end
