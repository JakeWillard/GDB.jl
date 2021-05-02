
abstract type Physical end

mutable struct Variable{T<:Physical}

    pval::Matrix{Float64} # past value
    cval::Matrix{Float64} # current value
    fval::Matrix{Float64} # future value
    x::Matrix{Float64}
    xx::Matrix{Float64}
    xy::Matrix{Float64}
    y::Matrix{Float64}
    yy::Matrix{Float64}
    t::Matrix{Float64}

end


function Variable{T}(pval::Matrix{Float64}, cval::Matrix{Float64}) where {T<:Physical}

    blank = zeros(Float64, size(pval))
    return Variable{T}(pval, cval, blank, blank, blank, blank, blank, blank, blank)
end



function Variable{T}(v0::Vector{Float64}, Nz::Int) where {T<:Physical}

    Nk = length(v0)
    pval = zeros(Float64, (Nk, Nz))
    cval = zeros(Float64, (Nk, Nz))

    for i=1:Nz
        pval[:,i] = v0[:]
        cval[:,i] = v0[:]
    end

    return Variable{T}(pval, cval)
end


function Variable{T}(dset) where {T<:Physical}

    if size(dset)[3] > 1
        return Variable{T}(dset[:,:,end-1], dset[:,:,end])
    elseif size(dset)[3] == 1
        return Variable{T}(dset[:,:,1], dset[:,:,1])
    end
end


function Variable{T}(f::Function, Nz::Int, grid::Grid) where {T<:Physical}

    Nk = grid.Nk
    v0 = zeros(Float64, Nk)

    for k=1:Nk
        x = grid.points[1,k]
        y = grid.points[2,k]
        v0[k] = f(x, y)
    end

    return Variable{T}(v0, Nz)
end


function timeshift!(var::Variable)

    var.pval[:,:] = var.cval[:,:]
    var.cval[:,:] = var.fval[:,:]
end


function to_3d_mesh(var::Variable, grid::Grid)

    vals = grid.inverse_projection * var.cval[:,:]

    Nx = grid.Nx
    Ny = grid.Ny
    Nz = size(vals)[2]
    out = zeros(Float64, (grid.Nx, grid.Ny, Nz))

    for j=1:Ny
        out[:,j,:] = vals[1+(j-1)*Nx:j*Nx,:]
    end

    return out
end