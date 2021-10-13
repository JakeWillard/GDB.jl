
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
        As = []
        for slc in picked
            append!(As, [sparse(view(A, :, slc))])
        end

        deltas = pmap(As) do M
            M \ r
        end

        for i=1:Np
            slc = picked[i]
            x[slc] = x[slc] + deltas[i]
        end

        r[:] = b - A*x
        err = norm(r) / bnorm
    end

    return x
end
