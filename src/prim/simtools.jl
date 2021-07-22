

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

function save_setup(path, stp::Setup)

    # open file
    fid = h5open(path, "w")

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

    fid["Grid/nanmask"] = grd.nanmask[:,:]
    fid["Grid/points"] = grd.points[:,:]
    fid["Grid/origin"] = grd.origin[:]
    fid["Grid/dims"] = Int64[grd.Nx, grd.Ny, grd.Nk, grd.mx, grd.my]
    fid["Grid/spacing"] = Float64[grd.dx, grd.dy]
    fid["Grid/Mxinv"] = grd.Mxinv[:,:]
    fid["Grid/Myinv"] = grd.Myinv[:,:]
    fid["Grid/MxyinvT"] = grd.MxyinvT[:,:]

    # save operators
    create_group(fid, "Operators")
    for (name, M) in stp.operators
        save_sparse_matrix(fid, M, "Operators/$(name)")
    end

    # close file
    close(fid)
end


function Setup(path)

    # open file
    fid = h5open(path, "r")

    # read grid
    projection = load_sparse_matrix(fid, "Grid/Proj")
    inverse_projection = load_sparse_matrix(fid, "Grid/ProjInv")
    restrictions = SparseMatrixCSC[]
    interpolations = SparseMatrixCSC[]

    Nl = length(keys(fid["Grid"])) - 7
    for i=1:Nl
        append!(restrictions, [load_sparse_matrix(fid, "Grid/Restr_$(i)")])
        append!(interpolations, [load_sparse_matrix(fid, "Grid/Interp_$(i)")])
    end

    nanmask = fid["Grid/nanmask"][:,:]
    points = fid["Grid/points"][:,:]
    origin = fid["Grid/origin"][:]
    Nx, Ny, Nk, mx, my = fid["Grid/dims"][:]
    dx, dy = fid["Grid/spacing"][:]
    Mxinv = fid["Grid/Mxinv"][:,:]
    Myinv = fid["Grid/Myinv"][:,:]
    MxyinvT = fid["Grid/MxyinvT"][:,:]
    grid = Grid(Grid(points, origin, projection, inverse_projection, restrictions, interpolations,
                     nanmask, mx, my, Nx, Ny, Nk, Mxinv, Myinv, MxyinvT, dx, dy))


    # read operators
    operators = Dict{String, SparseMatrixCSC}
    for name in keys(fid["Operators"])
        operators[name] = load_sparse_matrix(fid, "Operators/$(name)")
    end

    # close file
    close(fid)

    return Setup(grid, operators)
end


function save_variables(path::String, vars...)

    # open file
    fid = h5open(path, "r+")

    # read time index
    t = fid["t"]

    # save each var
    for var in vars
        fid["Output/$(var.name)"][:,:,t] = var[:,:,2]
    end

    # increase time index
    fid["t"] = t + 1

    # close file
    close(fid)

end


# NOTE: lets assume always that the first element of params is Nksp, the number of timesteps between saves.
function integrate_then_save(path::String, tstep::Function, params::Vector{Float64}, stp::Setup, vars...)

    Nskp = Int64(params[1])
    for _=1:Nskp
        tstep(params, stp, vars...)
    end

    save_variables(path, vars...)
end
