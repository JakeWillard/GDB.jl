


# struct RandomPartitionedQR
#
#     M :: SparseMatrixCSC
#     factorizations :: Array{SuiteSparse.SPQR.QRSparse{Float64, Int64}}
#     slices :: Vector{UnitRange{Int64}}
#     Np :: Int64
#
# end
#
#
# function RandomPartitionedQR(A, Np; Nc=50)
#
#     Nk = size(A)[2]
#     cs1 = div(Nk, Nc)
#     cs2 = div(Nk, Nc+5)
#     cs3 = div(Nk, Nc-5)
#     slices = [collect(Iterators.partition(1:Nk, cs1))...,
#               collect(Iterators.partition(1:Nk, cs2))...,
#               collect(Iterators.partition(1:Nk, cs3))...]
#
#     factorizations = [qr(A[:,s]) for s in slices]
#     return RandomPartitionedQR(A, factorizations, slices, Np)
# end
#
#
# function Base.:\(A::RandomPartitionedQR, b::Vector{Float64})
#
#     Nk = length(b)
#     Nm = length(A.factorizations)
#
#     x = zeros(Nk)
#     r = b - A*x
#     bnorm = norm(b)
#     err = norm(r) / norm(b)
#
#     while err > 1e-10
#
#         # only check for convergence after every 10th iteration
#         for _=1:10
#             is = shuffle(1:Nm)[1:A.Np]
#             d = pmap(is) do i
#                 A.factorizations[i] \ r
#             end
#             for j=1:A.Np
#                 slc = A.slices[is[j]]
#                 x[slc] = x[slc] + d[j]
#             end
#             r[:] = b - A*x
#         end
#
#         err = norm(r) / bnorm
#     end
#
#     return x
# end

# parallel-random-block-gauss-seidel
function prbgs(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64}; block_fraction=0.2, Np=2)

    Nk = length(x0)
    chunksize = Int64(floor(block_fraction*Nk/Np))
    offset = div(chunksize, 2)
    slices = [collect(Iterators.partition(1:Nk, chunksize))...,
              collect(Iterators.partition(offset:Nk, chunksize))...]
    Ns = length(slices)

    x = x0[:]
    r = b - A*x
    bnorm = norm(b)

    err = norm(r) / bnorm
    while err > 1e-10

        picked = [slices[rand(1:end)] for _=1:Np]

        deltas = pmap(picked) do slc
            A[:,slc] \ r
        end

        for i=1:Np
            slc = picked[i]
            x[slc] = x[slc] + deltas[i]
        end

        r[:] = b - A*x
        err = norm(r) / bnorm
        @info err
    end

    return x
end
