
function local_relaxation(x::Vector{Float64}, M::SparseMatrixCSC, f::Vector{Float64}, i::Int32, n::Int32, chunksize::Int32)
    for dummy=1:n
        x[i:i+chunksize] = M[i:i+chunksize,:] * x + f[i:i+chunksize]
    end
    return x[i:i+chunksize]
end


function prelaxation!(x::Vector{Float64}, M::SparseMatrixCSC, f::Vector{Float64}, n::Int32, chunksize::Int32)

    x = vcat(pmap(x -> local_relaxation(x, M, f, i, n, chunksize), i=1:chunksize:size(x)[1]-chunksize))
end


function jacobi(A::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64}; a=1.0, N=100, n=10, chunks=5)

    Dinv = a*inv(Diagonal(A))
    L = -tril(A, -1)
    U = -triu(A, 1)

    M = Dinv * (L + U)
    f = Dinv * b

    # XXX: assuming for now that size(A) evenly divides chunks
    chunksize = Int(ceil(size(A)[1] / chunks))

    for dummy=1:N
        prelaxation!(x, M, f, n, chunksize)
    end

end
