

abstract type Vorticity <: Physical end
abstract type MagneticPotential <: Physical end
abstract type ElectricPotential <: Physical end
abstract type ElectronTemp <: Physical end


function make_grid(a::Float64, N::Int, m::Int)

    inside(x, y) = (sqrt((x-0.5)^2 + (y-0.5)^2) < a)

    return Grid(inside, N, N, m, m)

end


function build_matrices(grid::Grid)

    Dx = x_derivative(1, grid)
    Dxx = x_dervative(2, grid)
    Dy = y_derivative(1, grid)
    Dyy = y_derivative(2, grid)
    Dxy = Dx * Dy
    L = laplacian(grid)

    return SparseMatrixCSC[Dx, Dxx, Dy, Dyy, Dxy, L]
end


function Circle(a::Float64, l::Float64, grid::Grid)

    pen = zeros(Float64, grid.Nk)

    for k=1:grid.Nk
        x = grid.points[1,k]
        y = grid.points[2,k]
        r = sqrt(x^2 + y^2)
        pen[k] = 0.5*(1 + tanh((r-a)/l))
    end

    PEN = Diagonal(pen)



end
