
#NOTE: this code is used in the example programs found in the /models folder, but I don't think its a good way
#      to structure the programs anymore, it seems to overcomplicate things. 

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



function Variable{T}(dset::Array{Float64, 3}) where {T<:Physical}

    pval = dset[:,:,end-1]
    cval = dset[:,:,end]
    return Variable{T}(pval, cval)
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



function Variable{T}(f::Function, Nz::Int, grid::Grid) where {T<:Physical}

    v0 = f_to_grid(f, grid)
    return Variable{T}(v0, Nz)
end
