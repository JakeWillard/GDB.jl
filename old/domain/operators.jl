

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
    D2d = P * D * Pinv / ((grd.dx)^nx *(grd.dy)^ny)

    return kron(sparse(I, grd.Nz, grd.Nz), D2d)
end


function line_average(points::Matrix{Float64}, mx::Int64, my::Int64, MinvT::Matrix{Float64}, grd::Grid)

    N = size(points)[2]
    is = ones(Int64, N*mx*my)
    js = Int64[]
    dats = Float64[]

    for i=1:N
        row_dat, row_j = interpolation_row(points[:,i]..., mx, my, MinvT, grd)
        js = [js; row_j]
        dats = [dats; row_dat]
    end

    return sparse(is, js, dats, 1, grd._Nx*grd._Ny) * transpose(grd.Proj) / N
end
