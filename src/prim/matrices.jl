
function x_derivative(n, grid::Grid)

    Nx = grid.Nx
    Ny = grid.Ny
    N = Nx*Ny

    Dx = sparse(zeros(Float64, (N, N)))
    Iy = sparse(I, Ny, Ny)
    stencil = grid.x_stencil[n,:]
    m = length(stencil)
    ic = Int(ceil(m/2.0))

    for i=1:m
        k = i - ic
        vec = stencil[i] * ones(Nx - abs(k))
        Dx += spdiagm(k => vec)
    end

    D = kron(Iy, Dx)
    Pinv = grid.inverse_projection
    P = grid.projection

    return P * D * Pinv
end


function y_derivative(n, grid::Grid)

    Nx = grid.Nx
    Ny = grid.Ny
    N = Nx*Ny

    Dy = sparse(zeros(Float64, (N, N)))
    Ix = sparse(I, Nx, Nx)
    stencil = grid.y_stencil[n,:]
    m = length(stencil)
    ic = Int(ceil(m/2.0))

    for i=1:m
        k = i - ic
        vec = stencil[i] * ones(Ny - abs(k))
        Dy += spdiagm(k => vec)
    end

    D = kron(Dy, Ix)
    Pinv = grid.inverse_projection
    P = grid.projection

    return P * D * Pinv
end


function laplacian(grid::Grid)

    alpha = (grid.dx / grid.dy) ^2

    Dxx = x_derivative(2, grid)
    Dyy = y_derivative(2, grid)

    return Dxx + alpha * Dyy
end
