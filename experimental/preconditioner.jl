
# take matrix A and return a matrix with the same sparsity pattern as A, but with
# off-diagonal elements diminished by factor proportional to abs(i - j)
function bandwidth_preconditioner(A::SparseMatrixCSC, bandwidth::Int64)

    is, js, dat = findnz(A)
    n, m = size(A)
    nnz = length(dat)
    dat_pc = zeros(nnz)

    for k=1:nnz
        u = maximum([0, 1 - abs(is[k] - js[k]) / bandwidth])

    end

    return sparse(is, js, dat_pc, n, m)
end


function matrix_split(A::SparseMatrixCSC, delta::Int64)

    is, js, dat = findnz(A)
    n, m = size(A)
    nnz = length(dat)
    dat_M = zeros(nnz)
    dat_N = zeros(nnz)

    for k=1:nnz
        u = maximum([0, 1 - abs(is[k] - js[k]) / delta])
        dat_M[k] = dat[k] * u
        dat_N[k] = dat[k] * (u - 1)
    end

    M = sparse(is, js, dat_M, n, m)
    N = sparse(is, js, dat_N, n, m)
    return M, N
end


function optimal_splitting(A::SparseMatrixCSC)

    is, js, dat = findnz(A)
    n, m = size(A)
    rad(a) = begin
        dat_m = a .* dat
        dat_n = (a .- 1) .* dat
        M = sparse(is, js, a .* dat, n, m)
        N = sparse(is, js, (a .- 1) .* dat, n, m)
        C = inv(Matrix(M)) * Matrix(N)
        return opnorm(C)
    end

    a0 = ones(length(dat))
    return optimize(rad, a0)
end


function spec_rad(M::SparseMatrixCSC)

    n,_ = size(M)
    v = rand(n)
    v = v / norm(v)
    rho = 2
    for _=1:100
        u = M * v
        rho_next = dot(u, v)
        v = u / norm(u)

        if abs(rho_next - rho) < 1e-8
            rho = rho_next
            break
        end
        rho  = rho_next
    end

    return rho
end



function optimal_scaling(A::SparseMatrixCSC, t)

    n, _ = size(A) # assume square A
    adiag = diag(A)

    rad(a) = begin
        d1 = a[1:n]
        d2 = a[n+1:end]
        D = Diagonal(1 ./ (d1.*adiag.*d2))
        C = I - D*(Diagonal(d1)*A*Diagonal(d2))
        return spec_rad(C)
    end

    a0 = ones(2*n)
    res = optimize(rad, a0, NelderMead(), Optim.Options(time_limit=t))
    af = Optim.minimizer(res)
    d1 = af[1:n]
    d2 = af[n+1:end]

    D = Diagonal(1 ./ (d1.*adiag.*d2))
    C = I - D*(Diagonal(d1)*A*Diagonal(d2))
    Sl = Diagonal(d1)
    Sr = Diagonal(d2)
    Apc = Sl * A * Sr
    D = Diagonal(1 ./ (d1.*adiag.*d2))
    C = I - D*Apc

    return C, D, Sl, Sr, res
end





# this algorithm is apparently successful as a two-sided preconditioner according to https://arxiv.org/pdf/1610.03871.pdf
# and also discussed in https://arxiv.org/pdf/1110.2805.pdf
function sinkhorn_knopp(A::SparseMatrixCSC; err=1e-8)

    n, m = size(A)
    At = transpose(A)
    d1 = rand(n)
    d2 = rand(m)
    r1 = 1
    r2 = 1

    for _=1:10000
        d1 = 1 ./ (A * d2)
        d2 = 1 ./ (At * d1)
    end

    return Diagonal(d1) * A * Diagonal(d2)
end
