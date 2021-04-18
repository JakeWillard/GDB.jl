
# penalize vector in favor of another vector
function penalize(v1::Vector{Float64}, v2::Vector{Float64}, pen::Vector{Float64})

    return v1 .* (1 .- pen) .+ pen .* v2
end


# penalize operator in favor of another operator
function penalize(A::SparseMatrixCSC, B::SparseMatrixCSC, pen::Vector{Float64})

    P = Diagonal(pen)
    return (I - P) * A + P * B
end
