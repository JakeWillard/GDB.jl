

function lsqr_next!(x::Vector{Float64}, r::Vector{Float64}, A::SparseMatrixCSC, pslices)

    delta = pmap(pslices) do inds
        A[:,inds] \ r
    end

    x[:] = x[:] + vcat(delta...)
end


function lsqr_solve(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64}, chunks::Int64)

    Nk = length(x0)
    chunksize = div(Nk, chunks)
    pslices = collect(Iterators.partition(1:Nk, chunksize))

    bnorm = norm(b)
    r = b - A*x0
    x = x0[:]
    err = norm(r) / bnorm

    while err > 1e-8
        lsqr_next!(x, r, A, pslices)
        r[:] = b - A*x
        err = norm(r) / bnorm
        @info err
    end

    return x
end


function lsqr_solve(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64}, chunks::Int64, dn::Int64)

    Nk = length(x0)
    chunksize = div(Nk, chunks)
    pslices = collect(Iterators.partition(1:Nk, chunksize))
    tups = [(s[1], s[end]+dn) for s in pslices[1:end-1]]
    append!(tups, [(pslices[end][1], pslices[end[2]])])

    bnorm = norm(b)
    r = b - A*x0
    x = x0[:]
    err = norm(r) / bnorm

    while err > 1e-8

        delta = pmap(tups) do t

        end
