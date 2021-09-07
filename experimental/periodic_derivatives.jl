
function stencil1d(m::Int)

    ic = Int(ceil(m/2.0))
    M = zeros((m, m))
    for i=1:m
        for j=1:m
            M[i, j] = (i - ic)^(j-1) / factorial(j-1)
        end
    end
    return inv(M)
end


function periodic_projection(N, m)

    js = [N-m+1:N...]
    append!(js, [1:N...], [1:m...])

    P = sparse([1:N...], [m+1:m+N...], ones(N), N, N+2*m)
    PT = sparse([1:N+2*m...], js, ones(N+2*m), N+2*m, N)

    return P, PT
end


function periodic_derivative_1d(n, m, N)

    Nb = N + 2*m
    D = sparse(zeros(Float64, (Nb, Nb)))
    stencil = stencil1d(m)[n+1,:]
    ic = Int(ceil(m/2.0))

    for i=1:m
        k = i - ic
        vec = stencil[i] * ones(Nb - abs(k))
        D += spdiagm(k => vec)
    end

    P, PT = periodic_projection(N, m)
    return P * D * PT
end


function periodic_derivative(nx, ny, mx, my, Nx, Ny, dx, dy)

    Dx = periodic_derivative_1d(nx, mx, Nx) / dx^nx
    Dy = periodic_derivative_1d(ny, my, Ny) / dy^ny
    return kron(Dy, Dx)
end
