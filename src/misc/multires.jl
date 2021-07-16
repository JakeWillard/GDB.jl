
mutable struct MultiResMatrix <: AbstractMatrix{Float64}

    mats :: Array{SparseMatrixCSC{Float64, Int32}, 1}
    level :: Int32

    MultiResMatrix(mats) = MultiResMatrix(mats, 1)
end


# minimal overloading for indexing purposes
Base.size(m::MultiResMatrix) = Base.size(m.mats[m.level])
Base.getindex(m::MultiResMatrix, inds...) = Base.getindex(m.mats[m.level], inds...)
Base.setindex!(m::MultiResMatrix, X, inds...) = Base.setindex!(m.mats[m.level], X, inds...)

# functions to go up and down courseness levels
function uplevel!(M::MultiResMatrix)

    M.level += 1
end

function dwnlevel!(M::MultiResMatrix)

    M.level -= 1
end
