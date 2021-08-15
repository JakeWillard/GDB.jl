

# NOTE: this is mostly just copied from wikipedia
function gmres_cycle(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64}, m::Int64; err_thresh = 1e-20)

    # compute initial residual, norm(b) and initial error
    r = b - A * x0
    bnorm = norm(b)
    err = norm(r) / bnorm

    n = length(x0)
    Q = zeros(Float64, (n,m+1))
    H = zeros(Float64, (m+2,m+1))
    sn = zeros(Float64, m)
    cn = zeros(Float64, m)
    e1 = zeros(Float64, m+1)
    e1[1] = 1.0

    beta = norm(r) * e1
    Q[:,1] = r / norm(r)

    # begin iteration
    nk = 0
    for k=1:m
        nk += 1

        # do arnoldi iteration
        q = A * Q[:,k]
        for i=1:k
            H[i,k] = transpose(q) * Q[:,i]
            q = q - H[i,k] * Q[:,i]
        end
        qnorm = norm(q)
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
        err = abs(beta[k+1]) / bnorm

        if err < err_thresh
            break
        end
    end

    # compute correction to x
    y = H[1:nk,1:nk] \ beta[1:nk]
    x = x0 + Q[:,1:nk] * y
    return x, err
end


function gmres_solve(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64}, m::Int64; maxiters=10000, err_thresh = 1e-8)

    x = x0[:]
    err = 1.0
    iteration = 0

    while (err > err_thresh)
        x, err = gmres_cycle(A, x, b, m, err_thresh=err_thresh)
        iteration += 1

        if iteration > maxiters
            @warn "Iteration limit reached, defaulting to direct solver."
            return A \ b
        end
    end

    return x
end
