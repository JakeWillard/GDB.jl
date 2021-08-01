

struct LinearLeftHandSide

    A :: SparseMatrixCSC
    M :: SparseMatrixCSC
    D :: SparseMatrixCSC

    LinearLeftHandSide(A::SparseMatrixCSC, w::Float64) = new(A, compute_jacobi_matrices(A, w)...)
end




function jacobi_preconditioned_gmres(L::LinearLeftHandSide, x0::Vector{Float64}, b::Vector{Float64}; n=100, m=10, err_thresh=1e-10)

    # precodition using jacobi_relaxation()
    x = jacobi_relaxation(n, L.A, L.M, L.D, x0, b)

    # finish with gmres
    x = gmres_solve(L.A, x, b, m; err_thresh=err_thresh)

    return x
end


function parallel_jacobi_preconditioned_gmres(L::LinearLeftHandSide, x0::Vector{Float64}, b::Vector{Float64}, pchunks::Int64, Nz::Int64; nsynch=100, nasynch=5, m=10, err_thresh=1e-10)

    # jacobi preconditioning
    x = distribute(x0; dist=pchunks*Nz)
    bdist = distribute(b; dist=pchunks*Nz)
    for _=1:nsynch
        x = asynch_jacobi_relaxation(nasynch, L.M, L.D, x, bdist)
    end

    # do gmres in parallel
    x = DArray((size(x)[1],), workers()[1:Nz], Nz) do inds

        xs = Array(x0[inds[1],1])
        bs = b[inds[1]]
        gmres_solve(L.A, xs, bs, m, err_thresh=err_thresh)
    end

    return Array(x)
end
