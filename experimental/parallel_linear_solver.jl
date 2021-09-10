
function rgbs(A, x0, b, Nchunk, Npar)

    Nk = length(x0)
    chunksize = div(Nk, Nchunk)
    slc = collect(Iterators.partition(1:Nk, chunksize))
    matlist = [qr(A[:,s]) for s in slc]

    x = x0[:]
    r = b - A*x
    for _=1:100
        parts = shuffle(1:length(matlist))[1:Npar]
        deltas = pmap(parts) do i
            matlist[i] \ r
        end

        for i=1:Npar
            rng = slc[parts[i]]
            x[rng] = x[rng] + deltas[i]
        end
        r[:] = b - A*x

        @info norm(r) / norm(b)
    end

    return x
end


function random_row_preconditioner!(A::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64}, Niter::Int64, Npar::Int64, slc::Vector{UnitRange{Int64}})

    factorized = [qr(A[:,s]) for s in slc]
    Nslc = length(slc)
    r = b - A*x

    for _=1:Niter
        partitions = shuffle(1:Nslc)[1:Npar]

        deltas = pmap(partitions) do i
            factorized[i] \ r
        end

        for j=1:Npar
            i = partitions[j]
            rng = slc[i]
            x[rng] = x[rng] + deltas[j]
        end

        r[:] = b - A*x
    end
end






function block_jacobi_preconditioner!(A::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64}, Nout::Int64, Nin::Int64, slices...)

    Nk = length(x)
    Nfac = length(slices)
    r = b - A*x
    factorized = [[qr(A[:,s]) for s in slc] for slc in slices]
    Adiag = diag(A)

    for _=1:Nout

        # d = zeros(Nk)
        # for qrlist in factorized
        #     d_slc = pmap(qrlist) do M
        #         M \ r
        #     end
        #     d += vcat(d_slc...) / Nfac
        # end
        #
        # x[:] = x + vcat(d)
        # r[:] = b - A*x

        # do pure jacobi smoothing
        for _=1:Nin
            x[:] = x + r ./ Adiag
            r[:] = b - A*x
        end
    end
end


function pmatmul(A::SparseMatrixCSC, v::Vector{Float64}, slc::Vector{UnitRange{Int64}})

    x = @distributed (+) for s in slc
        A[:,s] * v[s]
    end

    return x
end


function pdot(u::Vector{Float64}, v::Vector{Float64}, slc::Vector{UnitRange{Int64}})

    x = @distributed (+) for s in slc
        dot(u[s], v[s])
    end

    return x
end


function pnorm(u::Vector{Float64}, slc::Vector{UnitRange{Int64}})

    return sqrt(pdot(u, u, slc))
end


function pgmres(A::SparseMatrixCSC, x::Vector{Float64}, r::Vector{Float64}, m::Int64, slc::Vector{UnitRange{Int64}})

    n = length(x0)
    Q = zeros(Float64, (n,m+1))
    H = zeros(Float64, (m+2,m+1))
    sn = zeros(Float64, m)
    cn = zeros(Float64, m)
    e1 = zeros(Float64, m+1)
    e1[1] = 1.0

    rnorm = pnorm(r, slc)
    beta = rnorm* e1
    Q[:,1] = r / rnorm

    # begin iteration
    nk = 0
    for k=1:m
        nk += 1

        # do arnoldi iteration
        q = pmatmul(A, Q[:,k], slc)
        for i=1:k
            H[i,k] = pdot(q, Q[:,i], slc)
            q = q - H[i,k] * Q[:,i]
        end
        qnorm = pnorm(q, slc)
        H[k+1,k] = qnorm
        Q[:,k+1] = q / qnorm

        # apply givens rotation
        for i=1:k-1
            temp = cn[i] * H[i,k] + sn[i] * H[i+1,k]
            H[i+1,k] = -sn[i] * H[i,k] + cn[i] * H[i+1,k]
            H[i,k] = temp
        end

        # update sn and cn
        t = sqrt(H[k,k]^2 + H[k+1,k]^2)
        cn[k] = H[k,k] / t
        sn[k] = H[k+1,k] / t

        # eliminate H[k+1, k]
        H[k,k] = cn[k] * H[k,k] + sn[k] * H[k+1,k]
        H[k+1,k] = 0.0

        # update residual
        beta[k+1] = -sn[k]*beta[k]
        beta[k] = cn[k]*beta[k]
    end

    # compute correction to x
    y = H[1:nk,1:nk] \ beta[1:nk]
    return x + Q[:,1:nk] * y
end


function plinsolve(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64}, Nc::Int64; Np=10, Npar=5, Ng=10)

    Nk = length(x0)
    chunksize = div(Nk, Nc)
    slc = collect(Iterators.partition(1:Nk, chunksize))

    x = x0[:]
    random_row_preconditioner!(A, x, b, Np, Npar, slc)

    r = b - A*x
    # return r
    bnorm = norm(b)
    err = norm(r) / bnorm
    @info err

    while err > 1e-8
        x[:] = pgmres(A, x, r, Ng, slc)
        r[:] = b - A*x
        err = norm(r) / bnorm
        @info err
    end

    return x
end
