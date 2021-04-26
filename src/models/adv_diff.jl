

abstract type StreamFunction <: Physical end
abstract type Density <: Physical end


function BasicGrid(N::Int, m::Int)

    inside(x, y) = true

    return Grid(inside, N, N, m, m)

end



function SquareEdges(Ng, dn, N)

    R = sparse(zeros(Int32, (N, N)))
    for i=1:Ng
        k = 2*Ng - i
        R[i,k] = 1
    end
    for i=N-Ng:N
        k=2*(N-Ng) - i
        R[i,k] = 1
    end

    REF = kron(R, R)


    pen = zeros(N*N)
    for i=1:N
        for j=1:N
            k = i + N*(j-1)
            pen1 = 0.5*(1 + tanh((Ng-i)/dn))
            pen2 = 0.5*(1 + tanh((i - (N-Ng))/dn))
            pen3 = 0.5*(1 + tanh((Ng-j)/dn))
            pen4 = 0.5*(1 + tanh((j - (N-Ng))/dn))
            pen[k] = pen1 * pen2 * pen3 * pen4
        end
    end

    return Boundary(pen, Diagonal(pen), REF)
end


function partial_x!(f::Variable, t, Dx, dx)

    f.x = Dx * f.value[:,:,t] / dx
end


function partial_y!(f::Variable, t, Dy, dy)

    f.y = Dy * f.value[:,:,t] / dy
end



function partial_t!(n::Variable{Density}, phi::Variable{StreamFunction}, t, b::Boundary)

    n_t = phi.x .* n.y .- phi.y .* n.x
    n_ref = b.REF * n.value[:,:,t]
    n_t = penalize(n_t, n.value[:,:,t] .+ n_ref, b)
    n.t = n_t
end


function apply_diffusion!(n::Variable{Density}, L, dt, D, bdry::Boundary)

    n0 = n.value[:,1,3]
    A = penalize(I - dt*D*L, bdry.REF + I, bdry)
    b = penalize(n0, zeros(Float64, length(n0)), bdry.pen)

    n.value[:,1,3] = p_jacobi(A, n0, b, 0.2, 100, 1, 2, 1e-8)
end


function timestep!(n, phi, Dx, Dy, L, dt, D, bdry, dx, dy)

    # predictor step
    partial_x!(n, 2, Dx, dx)
    partial_y!(n, 2, Dy, dy)
    partial_t!(n, phi, 2, bdry)
    n.value[:,:,3] = 0.5*(n.value[:,:,1] .+ n.value[:,:,2]) .+ dt .* n.t
    apply_diffusion!(n, L, dt, D, bdry)

    # corrector step
    partial_x!(n, 3, Dx, dx)
    partial_y!(n, 3, Dy, dy)
    partial_t!(n, phi, 3, bdry)
    n.value[:,:,3] = n.value[:,:,2] .+ dt .* n.t
    apply_diffusion!(n, L, dt, D, bdry)

    n.value[:,:,1] = n.value[:,:,2]
    n.value[:,:,2] = n.value[:,:,3]
end


function simulation(n0::Function, phi0::Function, Nt, N, Ng)

    grid = BasicGrid(N, 5)
    n = Variable{Density}(n0, 1, grid)
    phi = Variable{StreamFunction}(phi0, 1, grid)

    Dx = x_derivative(2, grid)
    Dy = y_derivative(2, grid)
    L = laplacian(grid)

    partial_x!(phi, 2, Dx, grid.dx)
    partial_y!(phi, 2, Dy, grid.dy)

    wall = SquareEdges(Ng, Int(ceil(Ng/3)), N)

    dt = 0.01
    D = 0.5

    for dummy=1:Nt
        timestep!(n, phi, Dx, Dy, L, dt, D, wall, grid.dx, grid.dy)
    end

    return to_readable(n, grid)
end
