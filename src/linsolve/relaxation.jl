
function compute_jacobi_matrices(A::SparseMatrixCSC, w::Float64)

    # D = sparse(w*inv(Diagonal(A)))
    D = w*sparse(Diagonal(1 ./ diag(A)))
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



function asynch_jacobi_relaxation(N::Int64, M::SparseMatrixCSC, D::SparseMatrixCSC, x0::DArray, b::DArray)

    DArray((size(x0)[1],), procs(x0), length(procs(x0))) do inds

        Ms = M[inds[1],:]
        Ds = D[inds[1],:]
        xs = Array(x0)[:,1]
        fs = Ds * Array(b)[:,1]
        for _=1:N
            xs[inds[1]] = Ms * xs + fs
        end
        xs[inds[1],1]
    end

end
