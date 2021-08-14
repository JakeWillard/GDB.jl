


function jacobi_preconditioned_gmres(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64}, Nz, pchunks; w=2/3, Na=100, Ns=10, m=20, err_thresh=1e-8)

    Nk = div(length(x0), Nz)

    J = JacobiSmoother(A, b::Vector{Float64}, w, pchunks, Nk, Nz, Na)
    x = x0[:]
    x[:] = J*x
    for t=1:Ns
        x[:] = J * x
    end

    err = norm(A*x - b)
    if err > err_thresh
        zslices = collect(Iterators.partition(1:Nk*Nz, Nz))
        c = pmap(zslices) do inds
            Aloc = A[inds, inds]
            bloc = b[inds]
            gmres_solve(Aloc, x[inds], bloc, m; err_thresh=err_thresh)
        end
        return vcat(c...)
    else
        return x
    end

end





#
#
# struct LinearLeftHandSide
#
#     A :: SparseMatrixCSC
#     M :: SparseMatrixCSC
#     D :: SparseMatrixCSC
#
# end
# LinearLeftHandSide(A::SparseMatrixCSC, w::Float64) = LinearLeftHandSide(A, compute_jacobi_matrices(A, w)...)
#
#
#
# function jacobi_preconditioned_gmres(L::LinearLeftHandSide, x0::Vector{Float64}, b::Vector{Float64}; n=100, m=10, err_thresh=1e-10)
#
#     # precodition using jacobi_relaxation()
#     x = jacobi_relaxation(n, L.A, L.M, L.D, x0, b)
#
#     # finish with gmres
#     x = gmres_solve(L.A, x, b, m; err_thresh=err_thresh)
#
#     return x
# end
#
#
# function parallel_jacobi_preconditioned_gmres(L::LinearLeftHandSide, x0::Vector{Float64}, b::Vector{Float64}, pchunks::Int64, Nz::Int64; nsynch=100, nasynch=5, m=10, err_thresh=1e-10)
#
#     # jacobi preconditioning
#     x = distribute(x0; dist=pchunks*Nz)
#     bdist = distribute(b; dist=pchunks*Nz)
#     for _=1:nsynch
#         x = asynch_jacobi_relaxation(nasynch, L.M, L.D, x, bdist)
#     end
#
#     # do gmres in parallel
#     x = DArray((size(x)[1],), workers()[1:Nz], Nz) do inds
#
#         xs = Array(x0[inds[1],1])
#         bs = b[inds[1]]
#         gmres_solve(L.A, xs, bs, m, err_thresh=err_thresh)
#     end
#
#     return Array(x)
# end
