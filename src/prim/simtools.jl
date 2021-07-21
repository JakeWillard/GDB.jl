

abstract type Variable <: AbstractArray{Float64, 3} end

macro variable(sym)

    name = String(sym)
    type_def = Meta.parse("""
    mutable struct $(name) <: Variable
        vals :: Array{Float64, 3}
        name :: String

        Base.size(v::$(name)) = Base.size(v.vals)
        Base.getindex(v::$(name), inds...) = Base.getindex(v.vals, inds...)
        Base.setindex!(v::$(name), X, inds...) = Base.setindex!(v.vals, X, inds...)

        function $(name)(x::Matrix{Float64})
            n, m = size(x)
            X = zeros(Float64, (n, m, 2))
            X[:,:,1] = x[:,:]
            X[:,:,2] = x[:,:]
            new(X, "$(name)")
        end
    end
    """)

    return type_def
end


struct Setup

    grid::Grid
    operators::Dict{String, SparseMatrixCSC}

end


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
    grd = stp.grid
    create_group(fid, "Grid")
    save_sparse_matrix(fid, grd.projection, "Grid/Proj")
    save_sparse_matrix(fid, grd.inverse_projection, "Grid/ProjInv")

    for i=1:length(grd.restrictions)
        Restr = grd.restrictions[i]
        Interp = grd.interpolations[i]
        save_sparse_matrix(fid, Restr, "Grid/Restr_$(i)")
        save_sparse_matrix(fid, Restr, "Grid/Interp_$(i)")
    end

    fid["Grid/points"] = grd.points[:,:]
    fid["Grid/origin"] = grd.origin[:]
    fid["Grid/dims"] = Float64[grd.Nx, grd.Ny, grd.Nk, grd.mx, grd.my]
    fid["Grid/spacing"] = Float64[grd.dx, grd.dy]
    fid["Grid/Mxinv"] = grd.Mxinv[:,:]
    fid["Grid/Myinv"] = grd.Myinv[:,:]
    fid["Grid/MxyinvT"] = grd.MxyinvT[:,:]

    # save operators
    op = stp.operators
    create_group(fid, "Operators")
    for (name, M) in stp.operators
        save_sparse_matrix(fid, M, "Operators/$(name)")
    end
end
