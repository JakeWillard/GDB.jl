

function save_sparse_matrix(fid, A::SparseMatrixCSC, name::String)

    n, m = size(A)
    Is, Js, Vs = findnz(A)

    fid["$(name)_Is"] = Int32[Is...]
    fid["$(name)_Js"] = Int32[Js...]
    fid["$(name)_Vs"] = Float64[Vs...]
    fid["$(name)_size"] = Int32[n, m]
end


function load_sparse_matrix(fid, name::String)

    Is = fid["$(name)_Is"][:]
    Js = fid["$(name)_Js"][:]
    Vs = fid["$(name)_Vs"][:]
    n, m = fid["$(name)_size"][:]

    return sparse(Is, Js, Vs, n, m)
end


function save_grid(fid, name, grd::Grid)

    # create group for grid
    create_group(fid, name)

    # save projection
    save_sparse_matrix(fid, grd.Proj, "$(name)/Proj")
    fid["$(name)/r0"] = grd.r0[:]
    fid["$(name)/points"] = grd.points[:,:]
    fid["$(name)/spacing"] = Float64[grd.dx, grd.dy]
    fid["$(name)/_size"] = Int64[grd._Nx, grd._Ny]
    fid["$(name)/_NaNs"] = grd._nan_outside_boundaries[:,:]

end


function load_grid(fid, name)

    Proj = load_sparse_matrix(fid, "$(name)/Proj")
    r0 = fid["$(name)/r0"][:]
    points = fid["$(name)/points"][:,:]
    dx, dy = fid["$(name)/spacing"][:]
    Nx, Ny = fid["$(name)/_size"][:]
    NaNs = fid["$(name)/_NaNs"][:,:]

    return Grid(r0, points, Proj, dx, dy, Nx, Ny, NaNs)
end
