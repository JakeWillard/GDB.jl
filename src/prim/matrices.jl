
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



function reflection_matrix(delta, wall::Wall, grd::Grid)

    points = zeros(Float64, (2, grd.Nk))

    for k=1:grd.Nk
        x, y = grd.points[:,k]
        if 0 < smoothstep(x, y, delta, wall) < 1
            points[:,k] = wall.reflect(x, y)
        else
            points[:,k] = grd.points[:,k]
        end
    end

    return interpolation_matrix(points, grd)
end


function penalization_matrix(delta, wall::Wall, grd::Grid)

    return Diagonal(f_to_grid((x,y) -> smoothstep(x,y,delta,wall), grd))
end


mutable struct MultiResMatrix <: AbstractMatrix{Float64}

    mats :: Array{SparseMatrixCSC{Float64, Int32}, 1}
    lvl :: Int32

    MultiResMatrix(mats) = MultiResMatrix(mats, 1)
    Base.size(m::MultiResMatrix) = Base.size(m.mats[m.level])
    Base.getindex(m::MultiResMatrix, inds...) = Base.getindex(m.mats[m.level], inds...)
    Base.setindex!(m::MultiResMatrix, X, inds...) = Base.setindex!(m.mats[m.level], X, inds...)
    uplvl!(m::MultiResMatrix) = begin m.lvl += 1 end
    dwnlvl!(m::MultiResMatrix) = begin m.lvl -= 1 end
end


function galerkin_projection(A, grd)

    As = [A]
    for l=1:length(grd.restrictions)
        Restr = grd.restrictions[l]
        Interp = grd.interpolations[end-(l-1)]
        append!(As, Restr * As[end] * Interp)
    end

    return MultiResMatrix(As)
end


function interp1d(N, m)

    ic = Int(ceil(m/2.0))
    Minv = stencil1d(m)
    sten = transpose(Minv) * Float64[(1.5 - ic)^(j-1)/factorial(j-1) for j=1:m]

    Ih = zeros(Float64, (N,N))
    Is = I + Ih
    for i=1:m
        k = i - ic
        vec = sten[i] * ones(N - abs(k))
        Ih += diagm(k => vec)
    end

    Iout = zeros(Float64, (2*N,N))
    Iout[1:2:end,:] = Is
    Iout[2:2:end,:] = Ih
    return Iout
end
