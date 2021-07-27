

# functions for derivative matrices
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


# NOTE: Minv is the output of running stencil1d()
function regular_derivative_1d(N, n, Minv)

    D = sparse(zeros(Float64, (N, N)))
    stencil = Minv[n+1,:]
    m = length(stencil)
    ic = Int(ceil(m/2.0))

    for i=1:m
        k = i - ic
        vec = stencil[i] * ones(N - abs(k))
        D += spdiagm(k => vec)
    end

    return D
end


function derivative_matrix(nx, ny, Mxinv, Myinv, grd::Grid)

    Dx = regular_derivative_1d(grd._Nx, nx, Mxinv)
    Dy = regular_derivative_1d(grd._Ny, ny, Myinv)

    D = kron(Dy, Dx)
    P = grd.Proj
    Pinv = transpose(grd.Proj)

    return P * D * Pinv / ((grd.dx)^nx *(grd.dy)^ny)
end
