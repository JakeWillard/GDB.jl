
# take matrix A and return a matrix with the same sparsity pattern as A, but with
# off-diagonal elements diminished by factor proportional to abs(i - j)
function compute_preconditioner(A::SparseMatrixCSC, bandwidth::Int64)

    is, js, dat = findnz(A)
    n, m = size(A)
    nnz = length(dat)
    dat_pc = zeros(nnz)

    for k=1:nnz
        u = 1 - abs(is[k] - js[k]) / bandwidth
        if u > 0
            dat_pc[k] = dat[k] * (3*u^2 - 2*u^3)
        end
    end

    return sparse(is, js, dat_pc, n, m)
end
