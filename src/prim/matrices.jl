
function x_derivative(n, grid::Grid)

    Nx = grid.Nx
    Ny = grid.Ny

    Dx = sparse(zeros(Float64, (Nx, Nx)))
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

    Dy = sparse(zeros(Float64, (Ny, Ny)))
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


function interpolation_row(x, y, grid::Grid)

    Minv = grid.xy_stencil
    mx = size(grid.x_stencil)[1]
    my = size(grid.y_steincil)[2]

    x_ind = div(x, grid.dx)
    y_ind = div(y, grid.dy)
    xr = (x - grid.dx*x_ind) / grid.dx
    yr = (y - grid.dy*y_ind) / grid.dy

    v = zeros(Float64, (1, mx*my))
    for i=1:mx
        for j=1:my
            a = xr^(i-1) / factorial(i-1)
            b = yr^(j-1) / factorial(j-1)
            v[1,i + mx*(j-1)] = a*b
        end
    end

    dat = v * Minv

end
