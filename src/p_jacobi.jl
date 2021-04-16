
function relax(M, x, f, i, j, n)

    out = zeros(Float64, size(x))
    for dummy=1:n-1
        x[i:j] = M[i:j,:] * x + f[i:j]
    end
    out[i:j] = M[i:j,:] * x + f[i:j]
    return sparse(out)
end


function dist_jacobi(A::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64}; a=0.2, N=100, n=100, chunks=5)

    Dinv = a*inv(Diagonal(A))
    L = -tril(A, -1)
    U = -triu(A, 1)

    M = Dinv * (L + U)
    f = Dinv * b

    nodes = Int32[ceil(i) for i in LinRange(1, length(x), chunks)]
    xnext = zeros(Float64, size(x))

    err = 1.0
    while err > 1e-8
        for dummy=1:N
            # x = @distributed (+) for i=1:chunks-1
            #     relax(M, x, f, nodes[i], nodes[i+1], n)
            # end
            xnext[:] = sum(pmap(i -> relax(M, x, f, nodes[i], nodes[i+1], n), 1:chunks-1))
        end
        err = norm(xnext - x)
        println(err)
        x[:] = xnext[:]
    end
    return x


end
