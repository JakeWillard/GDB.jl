

function save_sparse_matrix(fid, A, name)

    n, m = size(A)
    Is, Js, Vs = findnz(A)

    fid["$(name)_Is"] = Int32[Is...]
    fid["$(name)_Js"] = Int32[Js...]
    fid["$(name)_Vs"] = Float64[Vs...]
    fid["$(name)_size"] = Int32[n, m]
end


function load_sparse_matrix(fid, name)

    Is = fid["$(name)_Is"][:]
    Js = fid["$(name)_Js"][:]
    Vs = fid["$(name)_Vs"][:]
    n, m = fid["$(name)_size"][:]

    return sparse(Is, Js, Vs, n, m)
end


function save_setup(fid, stp::Setup)

    # save grid
    grd = stp.grd
    create_group(fid, "Grid")
    save_sparse_matrix(fid, grd.projection, "Grid/Proj")
    save_sparse_matrix(fid, grd.inverse_projection, "Grid/ProjInv")
    fid["Grid/points"] = grd.points[:,:]
    fid["Grid/origin"] = grd.origin[:]
    fid["Grid/dims"] = Float64[grd.Nx, grd.Ny, grd.Nk, grd.mx, grd.my]
    fid["Grid/spacing"] = Float64[grd.dx, grd.dy]
    fid["Grid/Mxinv"] = grd.Mxinv[:,:]
    fid["Grid/Myinv"] = grd.Myinv[:,:]
    fid["Grid/MxyinvT"] = grd.MxyinvT[:,:]

    # save other matrices
    save_sparse_matrix(fid, stp.P1, "P1")
    save_sparse_matrix(fid, stp.P2, "P2")
    save_sparse_matrix(fid, stp.P3, "P3")
    save_sparse_matrix(fid, stp.R1, "R1")
    save_sparse_matrix(fid, stp.R2, "R2")
    save_sparse_matrix(fid, stp.R3, "R3")
    save_sparse_matrix(fid, stp.Dx, "Dx")
    save_sparse_matrix(fid, stp.Dy, "Dy")
    save_sparse_matrix(fid, stp.L, "L")
end


function load_setup(fid)

    # load grid
    points = fid["Grid/points"][:,:]
    origin = fid["Grid/origin"][:]
    P = load_sparse_matrix(fid, "Grid/Proj")
    Pinv = load_sparse_matrix(fid, "Grid/ProjInv")
    Nx, Ny, Nk, mx, my = fid["Grid/dims"][:]
    dx, dy = fid["Grid/spacing"][:]
    Mxinv = fid["Grid/Mxinv"][:,:]
    Myinv = fid["Grid/Myinv"][:,:]
    MxyinvT = fid["Grid/MxyinvT"][:,:]
    grd = Grid(points, origin, P, Pinv, mx, my, Nx, Ny, Nk, Mxinv, Myinv, MxyinvT, dx, dy)

    # load matrices
    P1 = load_sparse_matrix(fid, "P1")
    P2 = load_sparse_matrix(fid, "P2")
    P3 = load_sparse_matrix(fid, "P3")
    R1 = load_sparse_matrix(fid, "R1")
    R2 = load_sparse_matrix(fid, "R2")
    R3 = load_sparse_matrix(fid, "R3")
    Dx = load_sparse_matrix(fid, "Dx")
    Dy = load_sparse_matrix(fid, "Dy")
    L = load_sparse_matrix(fid, "L")

    return Setup(grd, P1, P2, P3, R1, R2, R3, Dx, Dy, L)
end
