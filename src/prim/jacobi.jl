
function relax(M, x, f, i, j, n)

    out = zeros(Float64, size(x))
    y = zeros(Float64, size(x))
    y[:] = x[:]
    for dummy=1:n-1
        y[i:j] = M[i:j,:] * y + f[i:j]
    end
    out[i:j] = M[i:j,:] * y + f[i:j]
    return sparse(out)
end


function p_jacobi(A::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64}, w::Float64, N::Int, n::Int, chunks::Int, threshold::Float64)

    Dinv = inv(Diagonal(A))
    # L = -tril(A, -1)
    # U = -triu(A, 1)
    bnorm = norm(b)

    M = I - w*Dinv*A
    # M = Dinv * (L + U)
    f = w*Dinv * b

    nodes = Int32[ceil(i) for i in LinRange(1, length(x), chunks)]
    xnext = zeros(Float64, size(x))

    err = 1.0
    while err > threshold
        for dummy=1:N
            xnext[:] = sum(pmap(i -> relax(M, x, f, nodes[i], nodes[i+1], n), 1:chunks-1))
        end
        err = norm(A * xnext - b) / bnorm

        println(err)
        x[:] = xnext[:]
    end
    return x


end
