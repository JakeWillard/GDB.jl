

function stencil(m)

    ic = Int(floor(m/2.0))
    M = zeros((m, m))
    for i=1:m
        for j=1:m
            M[i, j] = (i - ic)^(j-1) / factorial(j-1)
        end
    end
    return inv(M)
end


function laplacian2d(N, m)

    stens = stencil(m)
    ic = Int(floor((m+1)/2.0))

    D1d = zeros(Float64, (N, N))
    D1d[1,1] = 1.0
    D1d[N,N] = 1.0
    D1d_sparse = sparse(D1d)

    for i=1:m
        k = i-ic
        vec = stens[3,i] * ones(N - abs(k) - 2)
        D1d_sparse[2:N-1,2:N-1] += spdiagm(k => vec)
    end

    Isp = sparse(I, N, N)
    return kron(Isp, D1d_sparse) + kron(D1d_sparse, Isp)
end
