



struct Boundary

    pen::Vector{Float64}
    PEN::Diagonal
    REF::SparseMatrixCSC

end


function Boundary(p::Function, r::Function, grid::Grid)

end


function penalize(v1::Vector{Float64}, v2::Vector{Float64}, pen::Vector{Float64})

    P = Diagonal(pen)
    return (I - P) * v1 + P * v2
end


function penalize(v1::Matrix{Float64}, v2::Matrix{Float64}, pen::Vector{Float64})

    P = Diagonal(pen)
    return (I - P) * v1 + P * v2
end


function penalize(A::SparseMatrixCSC, B::SparseMatrixCSC, pen::Vector{Float64})

    P = Diagonal(pen)
    return (I - P) * A + P * B
end


# penalization given boundary
function penalize(v1::Matrix{Float64}, v2::Matrix{Float64}, b::Boundary)

    return penalize(v1, v2, b.pen)
end


# penalization of operator given boundary
function penalize(A::SparseMatrixCSC, B::SparseMatrixCSC, b::Boundary)

    return (I - b.PEN) * A + b.PEN * B
end
