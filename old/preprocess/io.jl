

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
    fid["$(name)/points"] = grd.points[:,:]
    fid["$(name)/spacing"] = Float64[grd.dx, grd.dy]
    fid["$(name)/_size"] = Int64[grd.Nk, grd.Nz, grd._Nx, grd._Ny]
    fid["$(name)/_NaNs"] = grd._nan_outside_boundaries[:,:]

    @info "Grid '$(name)' saved to disk."

end


function load_grid(fid, name::String)

    Proj = load_sparse_matrix(fid, "$(name)/Proj")
    r0 = fid["$(name)/r0"][:]
    points = fid["$(name)/points"][:,:]
    dx, dy = fid["$(name)/spacing"][:]
    Nk, Nz, Nx, Ny = fid["$(name)/_size"][:]
    NaNs = fid["$(name)/_NaNs"][:,:]

    @info "Loaded Grid '$(name)'."

    return Grid(r0, points, Proj, dx, dy, Nk, Nz, Nx, Ny, NaNs)
end


function save_lhs(fid, L::LinearLeftHandSide, name::String)

    # create group
    create_group(fid, name)
    save_sparse_matrix(fid, L.A, "$(name)/A")
    save_sparse_matrix(fid, L.M, "$(name)/M")
    save_sparse_matrix(fid, L.D, "$(name)/D")

    @info "LHS '$(name)' saved to disk."

end


function load_lhs(fid, name::String)

    A = load_sparse_matrix(fid, "$(name)/A")
    M = load_sparse_matrix(fid, "$(name)/M")
    D = load_sparse_matrix(fid, "$(name)/D")

    @info "Loaded LHS '$(name)'."
    return LinearLeftHandSide(A, M, D)
end


function save_assets(fid, a::Assets, name::String)

    # create group for assets
    create_group(fid, name)

    # save grid
    save_grid(fid, a.GRID, "$(name)/GRID")

    # save derivatives
    save_sparse_matrix(fid, a.Dx, "$(name)/Dx")
    save_sparse_matrix(fid, a.Dy, "$(name)/Dy")
    save_sparse_matrix(fid, a.Dxy, "$(name)/Dxy")
    save_sparse_matrix(fid, a.Dxx, "$(name)/Dxx")
    save_sparse_matrix(fid, a.Dyy, "$(name)/Dyy")
    save_sparse_matrix(fid, a.Dxxx, "$(name)/Dxxx")
    save_sparse_matrix(fid, a.Dyyy, "$(name)/Dyyy")
    save_sparse_matrix(fid, a.Dxxy, "$(name)/Dxxy")
    save_sparse_matrix(fid, a.Dxyy, "$(name)/Dxyy")
    save_sparse_matrix(fid, a.Ds, "$(name)/Ds")
    save_sparse_matrix(fid, a.Dss, "$(name)/Dss")

    # save left-hand-sides
    save_lhs(fid, a.DIFF_lnn, "$(name)/DIFF_lnn")
    save_lhs(fid, a.DIFF_lnTe, "$(name)/DIFF_lnTe")
    save_lhs(fid, a.DIFF_lnTi, "$(name)/DIFF_lnTi")
    save_lhs(fid, a.DIFF_u, "$(name)/DIFF_u")
    save_lhs(fid, a.DIFF_w, "$(name)/DIFF_w")
    save_lhs(fid, a.DIFF_A, "$(name)/DIFF_A")
    save_lhs(fid, a.HHOLTZ, "$(name)/HHOLTZ")

    # save boundary operators
    save_sparse_matrix(fid, a.P0, "$(name)/P0")
    save_sparse_matrix(fid, a.P1, "$(name)/P1")
    save_sparse_matrix(fid, a.P2, "$(name)/P2")
    save_sparse_matrix(fid, a.P3, "$(name)/P3")
    save_sparse_matrix(fid, a.R1, "$(name)/R1")
    save_sparse_matrix(fid, a.R2, "$(name)/R2")
    save_sparse_matrix(fid, a.R3, "$(name)/R3")
    save_sparse_matrix(fid, a.LAM, "$(name)/LAM")
    save_sparse_matrix(fid, a.DCHLT1, "$(name)/DCHLT1")
    save_sparse_matrix(fid, a.DCHLT2, "$(name)/DCHLT2")
    save_sparse_matrix(fid, a.DCHLT3, "$(name)/DCHLT3")
    save_sparse_matrix(fid, a.NMANN1, "$(name)/NMANN1")
    save_sparse_matrix(fid, a.NMANN2, "$(name)/NMANN2")
    save_sparse_matrix(fid, a.NMANN3, "$(name)/NMANN3")

    # save remaining assorted data
    save_sparse_matrix(fid, a.FLXAVG, "($name)/FLXAVG")
    fid["$(name)/TRGT"] = a.TRGT[:]
    fid["$(name)/Sn"] = a.Sn[:]
    fid["$(name)/STe"] = a.STe[:]
    fid["$(name)/STi"] = a.STi[:]
    fid["$(name)/params"] = a.params[:]
    fid["$(name)/dt"] = Float64[a.dt]
    fid["$(name)/N_subcycle"] = Int64[a.N_subcycle]
    # fid["$(name)/seed"][:] = a.seed[:]

    @info "Assets '$(name)' saved to disk."
end


function load_assets(fid, name::String)

    # load grid
    grd = load_grid(fid, "$(name)/GRID")

    # load derivatives
    Dx = load_sparse_matrix(fid, "$(name)/Dx")
    Dy = load_sparse_matrix(fid, "$(name)/Dy")
    Dxy = load_sparse_matrix(fid, "$(name)/Dxy")
    Dxx = load_sparse_matrix(fid, "$(name)/Dxx")
    Dyy = load_sparse_matrix(fid, "$(name)/Dyy")
    Dxxx = load_sparse_matrix(fid, "$(name)/Dxxx")
    Dyyy = load_sparse_matrix(fid, "$(name)/Dyyy")
    Dxxy = load_sparse_matrix(fid, "$(name)/Dxxy")
    Dxyy = load_sparse_matrix(fid, "$(name)/Dxyy")
    Ds = load_sparse_matrix(fid, "$(name)/Ds")
    Dss = load_sparse_matrix(fid, "$(name)/Dss")

    # load left-hand-sides
    DIFF_lnn = load_lhs(fid, "$(name)/DIFF_lnn")
    DIFF_lnTe = load_lhs(fid, "$(name)/DIFF_lnTe")
    DIFF_lnTi = load_lhs(fid, "$(name)/DIFF_lnTi")
    DIFF_u = load_lhs(fid, "$(name)/DIFF_u")
    DIFF_w = load_lhs(fid, "$(name)/DIFF_w")
    DIFF_A = load_lhs(fid, "$(name)/DIFF_A")
    HHOLTZ = load_lhs(fid, "$(name)/HHOLTZ")

    # load boundary operators
    P0 = load_sparse_matrix(fid, "$(name)/P0")
    P1 = load_sparse_matrix(fid, "$(name)/P1")
    P2 = load_sparse_matrix(fid, "$(name)/P2")
    P3 = load_sparse_matrix(fid, "$(name)/P3")
    R1 = load_sparse_matrix(fid, "$(name)/R1")
    R2 = load_sparse_matrix(fid, "$(name)/R2")
    R3 = load_sparse_matrix(fid, "$(name)/R3")
    LAM = load_sparse_matrix(fid, "$(name)/LAM")
    DCHLT1 = load_sparse_matrix(fid, "$(name)/DCHLT1")
    DCHLT2 = load_sparse_matrix(fid, "$(name)/DCHLT2")
    DCHLT3 = load_sparse_matrix(fid, "$(name)/DCHLT3")
    NMANN1 = load_sparse_matrix(fid, "$(name)/NMANN1")
    NMANN2 = load_sparse_matrix(fid, "$(name)/NMANN2")
    NMANN3 = load_sparse_matrix(fid, "$(name)/NMANN3")

    # load remaining assorted data
    FLXAVG = load_sparse_matrix(fid, "($name)/FLXAVG")
    TRGT = fid["$(name)/TRGT"][:]
    Sn = fid["$(name)/Sn"][:]
    STe = fid["$(name)/STe"][:]
    STi = fid["$(name)/STi"][:]
    params = fid["$(name)/params"][:]
    dt = fid["$(name)/dt"][1]
    N_subcycle = fid["$(name)/N_subcycle"][1]
    # seed = fid["$(name)/seed"][:]

    @info "Loaded Assets '$(name)'."

    return Assets(grd, Dx, Dy, Dxy, Dxx, Dyy, Dxxx, Dyyy, Dxxy, Dxyy, Ds, Dss,
                     DIFF_lnn, DIFF_lnTe, DIFF_lnTi, DIFF_u, DIFF_w, DIFF_A,
                     HHOLTZ, P0, P1, P2, P3, R1, R2, R3, LAM, DCHLT1, DCHLT2,
                     DCHLT3, NMANN1, NMANN2, NMANN3, FLXAVG, TRGT, Sn, STe,
                     STi, params, dt, N_subcycle)#, seed)
end
