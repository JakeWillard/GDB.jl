

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


# function lsqr_solve(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64}, chunks::Int64)
#
#     # get partition ranges
#     Nk = length(x0)
#     chunksize = div(Nk, chunks)
#     p_ranges = collect(Iterators.partition(1:Nk, chunksize))
#
#     for i=2:length(p_ranges)-1
#         sr = p_ranges[i-1][1]:p_ranges[i+1][2]
#         dr = length(p_ranges[i-1])+1:length(p_ranges[i-1])+length(p_ranges[i])
#         println(sr)
#     end
#     @assert false
#
#     # make list of ranges for subproblems
#     sub_ranges = [p_ranges[1][1]:p_ranges[3][2]]
#     append!(sub_ranges, [p_ranges[i-1][1]:p_ranges[i+1][2] for i=2:length(p_ranges)-1])
#     append!(sub_ranges, [p_ranges[end-2][1]:p_ranges[end][2]])
#
#     # make list of ranges for return deltas
#     d_ranges = [p_ranges[1]]
#     append!(d_ranges, [length(p_ranges[i-1])+1:(length(p_ranges[i-1])+length(p_ranges[i])) for i=2:length(p_ranges)-1])
#     append!(d_ranges, [(length(p_ranges[end-2])+length(p_ranges[end-1])):Nk])
#
#     bnorm = norm(b)
#     r = b - A*x0
#     x = x0[:]
#     err = norm(r) / bnorm
#
#     while err > 1e-8
#
#         delta = pmap(1:length(p_ranges)) do i
#             A[:,p_ranges[i]] \ r
#         end
#
#         x[:] = x + vcat([delta[i][d_ranges[i]] for i=1:length(p_ranges)]...)
#         r[:] = b - A*x
#         err = norm(r) / bnorm
#         @info err
#     end
#
#     return x
# end



# function lsqr_solve(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64}, chunks::Int64, dn::Int64)
#
#     Nk = length(x0)
#     chunksize = div(Nk, chunks)
#     pslices = collect(Iterators.partition(1:Nk, chunksize))
#     tups = [(s[1], s[end]+dn) for s in pslices[1:end-1]]
#     append!(tups, [(pslices[end][1], pslices[end[2]])])
#
#     bnorm = norm(b)
#     r = b - A*x0
#     x = x0[:]
#     err = norm(r) / bnorm
#
#     while err > 1e-8
#
#         delta = pmap(tups) do t
#
#         end
