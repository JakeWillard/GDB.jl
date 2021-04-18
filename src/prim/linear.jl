

function relax(M, x, f, i, j, n)

    out = zeros(Float64, size(x))
    for dummy=1:n-1
        x[i:j] = M[i:j,:] * x + f[i:j]
    end
    out[i:j] = M[i:j,:] * x + f[i:j]
    return sparse(out)
end


function p_jacobi(A::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64}, w::Float64, N::Int, n::Int, chunks::Int, threshold::Float64)

    Dinv = w*inv(Diagonal(A))
    L = -tril(A, -1)
    U = -triu(A, 1)

    M = Dinv * (L + U)
    f = Dinv * b

    nodes = Int32[ceil(i) for i in LinRange(1, length(x), chunks)]
    xnext = zeros(Float64, size(x))

    err = 1.0
    while err > threshold
        for dummy=1:N
            xnext[:] = sum(pmap(i -> relax(M, x, f, nodes[i], nodes[i+1], n), 1:chunks-1))
        end
        err = norm(xnext - x)
        x[:] = xnext[:]
        err = 0.0
    end
    return x


end


struct LinearSystem{T<:Physical}

    A::SparseMatrixCSC
    b::Vector{Float64}
    weight::Float64
    N::Int
    n::Int
    chunks::Int
    threshold::Float64

end


function solve(x0::Vector{Float64}, lin::LinearSystem)

    return p_jacobi(lin.A, x0, lin.b, lin.weight, lin.N, lin.n, lin.chunks, lin.threshold)
end
