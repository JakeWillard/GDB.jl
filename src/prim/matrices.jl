
function x_derivative(n, grid::Grid)

    Nx = grid.Nx
    Ny = grid.Ny

    Dx = sparse(zeros(Float64, (Nx, Nx)))
    Iy = sparse(I, Ny, Ny)
    stencil = grid.Mxinv[n+1,:]
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
    stencil = grid.Myinv[n+1,:]
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

    return interpolation_row(x, y, grid.MxyinvT, grid.Nx, grid.Ny, grid.mx, grid.my)
end



function interpolation_matrix(points, grid::Grid)

    dat = Float64[]
    is = Int32[]
    js = Int32[]

    Np = size(points)[2]
    for k=1:Np
        x = points[1,k]
        y = points[2,k]
        row_dat, row_js = interpolation_row(x, y, grid)
        row_is = k*ones(Int, grid.mx*grid.my)

        dat = [dat; row_dat]
        is = [is; row_is]
        js = [js; row_js]
    end

    return sparse(is, js, dat, Np, grid.Nx*grid.Ny) * grid.inverse_projection
end



function reflection_matrix(wall:Wall, grd::Grid)

    points = zeros(Float64, grd.Nk)
    n = 0

    for k=1:Nk
        x, y = grd.points[:,k]
        if 0 < wall.sstep(x,y) < 1
            n += 1
            points[:,n] = wall.reflect(x, y)
        else
            points[:,n] = grd.points[:,k]
        end
    end

    return interpolation_matrix(points, grd)
end
