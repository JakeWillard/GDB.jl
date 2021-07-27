

struct LinearLeftHandSide

    A :: SparseMatrixCSC
    M :: SparseMatrixCSC
    D :: SparseMatrixCSC

    LinearLeftHandSide(A::SparseMatrixCSC, w::Float64) = new(A, compute_jacobi_matrices(A, w)...)
end




function jacobi_preconditioned_gmres(L::LinearLeftHandSide, x0::Vector{Float64}, b::Vector{Float64}; n=100, m=10, err_thresh=1e-20)

    # compute M and D
    M, D = L.M, L.D

    # precodition using jacobi_relaxation()
    x = jacobi_relaxation(n, L.A, L.M, L.D, x0, b)

    # finish with gmres
    x = gmres_solve(L.A, x, b, m; err_thresh=err_thresh)

    return x
end
