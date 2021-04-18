
# penalize vector in favor of another vector
function penalize(v1::Vector{Float64}, v2::Vector{Float64}, pen::Vector{Float64})

    return v1 .* (1 .- pen) .+ pen .* v2
end


# penalize operator in favor of another operator
function penalize(A::SparseMatrixCSC, B::SparseMatrixCSC, P::Vector{Float64})

    return (I - P) * A + P * B
end


# penalize vector in favor of multiple vectors
function penalize(a::Vector{Float64}, b::Array{Vector{Float64}, 2})

    for i=1:size(b)[2]
        a = penalize(a, b[:,i]...)
    end
    return a
end


# penalize operator in favor of multiple operators
function penalize(A::SparseMatrixCSC, B::Array{SparseMatrixCSC, 2})

    for i=1:size(B)[2]
        A = penalize(A, B[:,i]...)
    end
    return A
end
