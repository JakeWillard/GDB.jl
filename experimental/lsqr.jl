

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

    facs = [qr(A[:,s]) for s in pslices]
    errors = zeros(1000)


    for i=1:1000

        delta = pmap(facs) do Aloc
            Aloc \ r
        end
        x[:] = x + vcat(delta...)
        r[:] = b - A*x

        errors[i] = norm(r) / bnorm
    end

    return errors
end


function lsqr_serial(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64}, N::Int64)

    D = (2/3)*Diagonal(1 ./ diag(A))
    C = (I - D*A)
    f = D*b

    Nk = length(x0)
    chunksize = div(Nk, N)
    pslices = collect(Iterators.partition(1:Nk, chunksize))

    As = [A[:,s] for s in pslices]
    qrs = [qr(hcat(As[end], As[1], As[2]))]
    drng = [length(pslices[end])+1:length(pslices[end])+length(pslices[1])]
    for i=2:length(As)-1
        append!(qrs, [qr(hcat(As[i-1], As[i], As[i+1]))])
        _a = length(pslices[i-1])
        _b = length(pslices[i])
        append!(drng, [_a+1:_a+_b])
    end
    append!(qrs, [qr(hcat(As[end-1], As[end], As[1]))])
    append!(drng, [length(pslices[end-1])+1:length(pslices[end-1])+length(pslices[end])])

    # @assert sum([length(d) for d in drng]) == Nk

    bnorm = norm(b)
    r = b - A*x0
    x = x0[:]

    Niter = 100
    errors = zeros(Niter)
    errors[1] = norm(r) / bnorm
    residuals = zeros(Nk, Niter)
    residuals[:,1] = r[:]
    for i=2:Niter

        delta = pmap(1:length(pslices)) do i
            M = qrs[i]
            rng = drng[i]
            d = M \ r
            d[rng]
        end

        x[:] = x + (2/3)*vcat(delta...)

        # do inner jacobi iteration
        for _=1:2
            x[:] = C*x + f
        end

        r[:] = b - A*x
        errors[i] = norm(r) / bnorm
        residuals[:,i] = r[:]

        @info errors[i]
    end

    return errors, residuals
end


function lsqr_avg(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64}, N::Int64, n::Int64)

    Nk = length(x0)
    chunksize_1 = div(Nk, N)
    chunksize_2 = div(Nk, N-n)
    chunksize_3 = div(Nk, N+n)
    pslices_1 = collect(Iterators.partition(1:Nk, chunksize_1))
    pslices_2 = collect(Iterators.partition(1:Nk, chunksize_2))
    pslices_3 = collect(Iterators.partition(1:Nk, chunksize_3))

    qrs_1 = [qr(A[:,s]) for s in pslices_1]
    qrs_2 = [qr(A[:,s]) for s in pslices_2]
    qrs_3 = [qr(A[:,s]) for s in pslices_3]

    bnorm = norm(b)
    r = b - A*x0
    x = x0[:]

    Niters = 100
    res = zeros(Nk, Niters)
    res[:,1] = r[:]

    for i=2:Niters

        d1 = pmap(qrs_1) do M
            M \ r
        end
        d2 = pmap(qrs_2) do M
            M \ r
        end
        d3 = pmap(qrs_3) do M
            M \ r
        end

        delta = (vcat(d1...) + vcat(d2...) + vcat(d3...)) / 3.0
        x[:] = x + delta
        r[:] = b - A*x
        res[:,i] = r[:]
        @info norm(r) / bnorm
    end

    return res
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
