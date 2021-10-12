
struct LHSMatrix

    M :: SparseMatrixCSC{Float64, Int64}
    Mdist :: DArray{Float64, 2, SparseMatrixCSC{Float64, Int64}}

end


function LHSMatrix(A::SparseMatrixCSC{Float64, Int64}; prols=workers())

    Mdist = DArray(size(A), workers(), (1,length(workers()))) do inds
        j = inds[2]
        A[:,j]
    end

    return LHSMatrix(A, Mdist)
end


function solve_local(M::DArray{Float64, 2, SparseMatrixCSC{Float64, Int64}}, r::Vector{Float64})

    delta = M[:L] \ r
    inds = localindices(M)[2]

    return delta, inds
end


function dlinsolve(L::LHSMatrix, x0::Vector{Float64}, b::Vector{Float64}, Np)

    n = length(x0)
    pids = procs(L.Mdist)
    nworkers = length(pids)
    chunks = collect(Iterators.partition(1:nworkers, Np))
    slices = fetch.([@spawnat pid localindices(L.Mdist)[2] for pid in pids])

    r = b - L.M*x0
    x = x0[:]

    bnorm = norm(b)
    err = norm(r) / bnorm
    while err > 1e-10
        for chunk in chunks

            res = Future[]
            for pid in pids[chunk]
                fut = @spawnat pid solve_local(L.Mdist, r)
                append!(res, [fut])
            end

            for f in res
                delta, inds = fetch(f)
                x[inds] += delta
            end
            
            r[:] = b - L.M*x
        end

        err = norm(r) / bnorm
        println(err)
    end

    return x
end
