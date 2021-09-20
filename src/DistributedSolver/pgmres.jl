

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
