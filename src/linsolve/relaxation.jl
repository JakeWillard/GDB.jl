
function compute_jacobi_matrices(A::SparseMatrixCSC, w::Float64)

    D = sparse(w*inv(Diagonal(A)))
    M = I - D*A
    return M, D
end


function jacobi_relaxation(N::Int64, A::SparseMatrixCSC, M::SparseMatrixCSC, D::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64})

    r = b - A*x0
    x = zeros(Float64, length(r))
    f = D * r

    for _=1:N
        x[:] = M * x + f
    end

    return x0 + x
end
