

function block_jacobi_preconditioner!(A::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64}, N::Int64, slices...)

    Nk = length(x)
    Nfac = length(slices)
    r = b - A*x
    factorized = [[qr(A[:,s]) for s in slc] for slc in slices]

    for _=1:N

        d = zeros(Nk)
        for qrlist in factorized
            d_slc = pmap(qrlist) do M
                M \ r
            end
            d += vcat(d_slc...) / Nfac
        end

        x[:] = x + vcat(d)
        r[:] = b - A*x

        # do some richardson iteration?
        for _=1:10
            x[:] = x + (2/3)*r
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


function plinsolve(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64}, Nc::Int64; Npc=10, Ng=10)

    Nk = length(x0)
    chunksize0 = div(Nk, Nc)
    chunksize1 = div(Nk, Nc-1)
    slc0 = collect(Iterators.partition(1:Nk, chunksize0))
    slc1 = collect(Iterators.partition(1:Nk, chunksize1))

    x = x0[:]
    block_jacobi_preconditioner!(A, x, b, Npc, slc0, slc0, slc1)

    r = b - A*x
    return r
    bnorm = norm(b)
    err = norm(r) / bnorm
    @info err

    while err > 1e-8
        # x[:] = pgmres(A, x, r, Ng, slc0)
        x[:] = gmres_cycle(A, x, b, Ng)
        r[:] = b - A*x
        err = norm(r) / bnorm
        @info err
    end

    return x
end
