



struct Boundary

    pen::Vector{Float64}
    PEN::Diagonal
    REF::SparseMatrixCSC

end


function Boundary(p::Function, r::Function, grid::Grid)

end




# penalize vector in favor of another vector
function penalize(v1::Vector{Float64}, v2::Vector{Float64}, pen::Vector{Float64})

    return v1 .* (1 .- pen) .+ pen .* v2
end


# penalize operator in favor of another operator
function penalize(A::SparseMatrixCSC, B::SparseMatrixCSC, pen::Vector{Float64})

    P = Diagonal(pen)
    return (I - P) * A + P * B
end


# penalization given boundary
function penalize(v1::Vector{Float64}, v2::Vector{Float64}, b::Boundary)

    return penalize(v1, v2, b.pen)
end


# penalization of operator given boundary
function penalize(A::SparseMatrixCSC, B::SparseMatrixCSC, b::Boundary)

    return (I - b.PEN) * A + b.PEN * B
end
