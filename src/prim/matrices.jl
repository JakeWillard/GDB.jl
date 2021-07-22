
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

    return P * D * Pinv / (grid.dx)^n
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

    return P * D * Pinv / (grid.dy)^n
end


function laplacian(grid::Grid)

    return x_derivative(2, grid) + y_derivative(2, grid)
end



function interpolation_row(x, y, grid::Grid)

    return interpolation_row(x, y, grid.origin[1], grid.origin[2], grid.dx, grid.dy, grid.MxyinvT, grid.Nx, grid.Ny, grid.mx, grid.my)
end


function reflection_matrix(delta::Float64, wall::Wall, grd::Grid)

    # define the function we want to map onto grd (see grids.jl for how Base.map has been overloaded.)
    function f(x0, y0, grd)
        if 0 < smoothstep(x0, y0, delta, wall) < 1
            x, y =  wall.reflect(x0, y0)
            row_dat, row_js = interpolation_row(x, y, grd)
        else
            # we want identity, but sometimes floor() messes this up so we should nudge x and y a small amount
            # XXX: check this later to see if it causes issues, these points shouldn't matter very much so I'm guessing it's fine.
            x = x0 + 0.00001*grd.dx
            y = y0 + 0.00001*grd.dy
            row_dat, row_js = interpolation_row(x, y, grd)
        end
        return hcat(row_dat, row_js)
    end

    M = grid_map(f, (grd.mx*grd.my,2), grd)
    dat = M[:,1]
    js = Int64[M[:,2]...]
    is = vcat([k*ones(Int64, grd.mx*grd.my) for k=1:grd.Nk]...)

    return sparse(is, js, dat, grd.Nk, grd.Nx*grd.Ny) * grd.inverse_projection
end


function penalization_matrix(delta, wall::Wall, grd::Grid)

    return Diagonal(f_to_grid((x,y) -> smoothstep(x,y,delta,wall), grd))
end


# NOTE: tracing parallel or anti-parallel to b can be achieved with positive or negative ds
function fieldline_map_matrix(bx::Function, by::Function, bz::Function, ds::Float64, Nz::Int64, grd::Grid)

    deltaPhi = 2*pi / Nz

    # define function we want to map onto grd
    function f(x0, y0, grd)

        x, y, deltaS = trace_fieldline(x0, y0, bx, by, bz, ds, deltaPhi)
        row_dat, row_js = interpolation_row(x, y, grd)
        dS_vec = zeros(Float64, grd.mx*grd.my)
        dS_vec[1] = deltaS

        return hcat(row_dat, row_js, dS_vec)
    end

    M = grid_map(f, (grd.mx*grd.my,3), grd)
    dS = M[1:grd.mx*grd.my:end,3]
    dat = M[:,1]
    js = Int64[M[:,2]...]
    is = vcat([k*ones(Int64, grd.mx*grd.my) for k=1:grd.Nk]...)

    matrix = sparse(is, js, dat, grd.Nk, grd.Nx*grd.Ny) * grd.inverse_projection
    return matrix, dS
end
