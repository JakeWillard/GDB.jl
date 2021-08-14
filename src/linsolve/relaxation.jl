
struct JacobiSmoother

    C :: SparseMatrixCSC
    f :: Vector{Float64}
    asynch_steps :: Int64
    slices :: Vector{UnitRange}

end


function JacobiSmoother(A::SparseMatrixCSC, b::Vector{Float64}, w::Float64, pchunks, Nk, Nz, Na)

    D = w*Diagonal(1 ./ diag(A))
    C = I - D*A
    f = D * b

    chunksize = div(Nk, pchunks)
    @info chunksize
    pslices = collect(Iterators.partition(1:Nk, chunksize))
    # pslices = UnitRange{Int64}[g[1]:g[end] for g in pgroups]
    slices = UnitRange{Int64}[]
    for z=1:Nz
        slices = [slices; [s .+ (z-1)*Nk for s in pslices]]
    end

    return JacobiSmoother(C, f, Na, pslices)
end


function Base.:*(J::JacobiSmoother, v::Vector{Float64})

    out = pmap(J.slices) do inds
        x = v[:]
        C_loc = J.C[inds, :]
        f_loc = J.f[inds]

        for t=1:J.asynch_steps
            x[inds] = C_loc * x + f_loc
        end
        x[inds]
    end

    return vcat(out...)
end
