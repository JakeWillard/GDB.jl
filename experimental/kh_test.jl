# kelvin helmholtz (wakatani with no alpha terms)
using LinearAlgebra, SparseArrays, ProgressMeter, ForwardDiff, RecipesBase, Plots, Distributed

include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/ImmersedMirrors/mirror.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/ImmersedMirrors/grid.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/ImmersedMirrors/extrapolator.jl")

# code taken from test_schemes.jl
function derivative_1d(N, n, m)

    # make stencil
    ic = Int(ceil(m/2.0))
    M = zeros((m, m))
    for i=1:m
        for j=1:m
            M[i, j] = (i - ic)^(j-1) / factorial(j-1)
        end
    end
    stencil = inv(M)[n+1,:]

    D = sparse(zeros(Float64, (N, N)))
    for i=1:m
        k = i - ic
        vec = stencil[i] * ones(N - abs(k))
        D += spdiagm(k => vec)
    end

    return D
end


function mixed_derivative(nx, ny, m, grd)

    operator_to_grid(grd) do
        Dx = derivative_1d(grd._Nx, nx, m)
        Dy = derivative_1d(grd._Ny, ny, m)
        kron(Dy, Dx) / grd.dr^(nx + ny)
    end
end


function Laplacian(m, grd)
    Dxx = mixed_derivative(2, 0, m, grd)
    Dyy = mixed_derivative(0, 2, m, grd)

    return Dxx + Dyy
end


function hyperdiffusion(n, mu, grd)

    operator_to_grid(grd) do
        Dx = derivative_1d(grd._Nx, 2*n, 2*n+1)
        Dy = derivative_1d(grd._Ny, 2*n, 2*n+1)
        Ix = sparse(I, grd._Nx, grd._Nx)
        Iy = sparse(I, grd._Ny, grd._Ny)

        I + (-1)^n*mu*(kron(Iy, Dx) + kron(Dy, Ix))
    end

end


function direct_solve(A, b, bval, extr)

    Anew, bnew = extr(A, b, bval)
    bnew = transpose(A)*b
    Anew = transpose(A)*A
    return extr(Anew \ bnew, bval)
end


function kelvinhelmholtz_update!(w, phi, n, Dx, Dy, w_ex, phi_ex, n_ex, L, Dh, dt)
    # goes through and solves each value twice - once using first two columns, then purely with second column
    bval = zeros(size(w)[1])

    phi_x = Dx*phi[:,2]
    phi_y = Dy*phi[:,2]
    n_x = Dx*n[:,2]
    n_y = Dy*n[:,2]
    w_x = Dx*w[:,2]
    w_y = Dy*w[:,2]

    w_t = -(phi_x.*w_y - phi_y.*w_x)
    # w_t = zeros(size(phi_x))
    n_t = -(phi_x.*n_y - phi_y.*n_x)
    # n_t = zeros(size(phi_x))

    w[:,3] = (w[:,1] + w[:,2])/2 + dt*w_t
    n[:,3] = (n[:,1] + n[:,2])/2 + dt*n_t

    w[:,3] = direct_solve(Dh, w[:,3], bval, w_ex)
    n[:,3] = direct_solve(Dh, n[:,3], bval, n_ex)
    phi[:,3] = direct_solve(L, w[:,3], bval, phi_ex)

    phi_x = Dx*phi[:,3]
    phi_y = Dy*phi[:,3]
    n_x = Dx*n[:,3]
    n_y = Dy*n[:,3]
    w_x = Dx*w[:,3]
    w_y = Dy*w[:,3]

    w_t = -(phi_x.*w_y - phi_y.*w_x)
    n_t = -(phi_x.*n_y - phi_y.*n_x)

    w[:,3] = w[:,2] + dt*w_t
    n[:,3] = n[:,2] + dt*n_t

    w[:,3] = direct_solve(Dh, w[:,3], bval, w_ex)
    n[:,3] = direct_solve(Dh, n[:,3], bval, n_ex)
    phi[:,3] = direct_solve(L, w[:,3], bval, phi_ex)

    w[:,1:2] = w[:,2:3]
    n[:,1:2] = n[:,2:3]
    phi[:,1:2] = phi[:,2:3]

    return w, n, phi
end


#=
Adding parallel derivative of current
Current density given by eq 19 in Ben's GDB paper
jparallel = c/4pi * perpendicular laplacian of psi
How to calculate flux?
- input function
=#


function kelvinhelmholtzcurrent_update!(w, phi, n, psi, j, Dx, Dy, w_ex, phi_ex, n_ex, psi_ex, j_ex, L, Dh, dt, eta)
    # goes through and solves each value twice - once using first two columns, then purely with second column
    bval = zeros(size(w)[1])

    phi_x = Dx*phi[:,2]
    phi_y = Dy*phi[:,2]
    n_x = Dx*n[:,2]
    n_y = Dy*n[:,2]
    w_x = Dx*w[:,2]
    w_y = Dy*w[:,2]
    psi_x = Dx*psi[:,2]
    psi_y = Dy*psi[:,2]
    j_x = Dx*j[:,2]
    j_y = Dy*j[:,2]

    w_t = -(phi_x.*w_y - phi_y.*w_x) + (psi_x.*j_y - psi_y.*j_x)
    n_t = -(phi_x.*n_y - phi_y.*n_x)
    psi_t = -(phi_x.*psi_y - phi_y.*psi_x) + eta*j[:,2]

    w[:,3] = (w[:,1] + w[:,2])/2 + dt*w_t
    n[:,3] = (n[:,1] + n[:,2])/2 + dt*n_t
    psi[:,3] = (psi[:,1] + psi[:,2])/2 + dt*psi_t

    w[:,3] = direct_solve(Dh, w[:,3], bval, w_ex)
    n[:,3] = direct_solve(Dh, n[:,3], bval, n_ex)
    psi[:,3] = direct_solve(Dh, psi[:,3], bval, psi_ex)
    phi[:,3] = direct_solve(L, w[:,3], bval, phi_ex)
    j[:,3] = L*psi[:,3]

    phi_x = Dx*phi[:,3]
    phi_y = Dy*phi[:,3]
    n_x = Dx*n[:,3]
    n_y = Dy*n[:,3]
    w_x = Dx*w[:,3]
    w_y = Dy*w[:,3]
    psi_x = Dx*psi[:,3]
    psi_y = Dy*psi[:,3]
    j_x = Dx*j[:,3]
    j_y = Dy*j[:,3]

    w_t = -(phi_x.*w_y - phi_y.*w_x) + (psi_x.*j_y - psi_y.*j_x)
    n_t = -(phi_x.*n_y - phi_y.*n_x)
    psi_t = -(phi_x.*psi_y - phi_y.*psi_x) + eta*j[:,3]

    w[:,3] = w[:,2] + dt*w_t
    n[:,3] = n[:,2] + dt*n_t
    psi[:,3] = psi[:,2] + dt*psi_t

    w[:,3] = direct_solve(Dh, w[:,3], bval, w_ex)
    n[:,3] = direct_solve(Dh, n[:,3], bval, n_ex)
    psi[:,3] = direct_solve(Dh, psi[:,3], bval, psi_ex)
    phi[:,3] = direct_solve(L, w[:,3], bval, phi_ex)
    j[:,3] = L*psi[:,3]

    w[:,1:2] = w[:,2:3]
    n[:,1:2] = n[:,2:3]
    phi[:,1:2] = phi[:,2:3]
    psi[:,1:2] = psi[:,2:3]
    j[:,1:2] = j[:,2:3]

    return w, n, phi, psi, j
end


# MAIN
m = 5
dt = .01
N = 900
mu = .3
eta = .001

thetas = LinRange(0, 2*pi, 17)[1:16]
c1 = zeros(2, 16)
c2 = zeros(2, 16)
for i=1:16
  c1[:,i] = [cos(thetas[i]), sin(thetas[i])]
  c2[:,i] = .5*[cos(thetas[i]), -sin(thetas[i])]
end
M = Mirror([c1, c2])

r0 = [-1.1, -1.1]
r1 = [1.1, 1.1]
h = .1
grid = Grid(M, h, r0, r1, 8)

n(x,y) = sqrt(x^2+y^2) - 1/3*x^3
n0 = function_to_grid(n, grid)
zros = zeros(grid.Nk)
n0 = hcat(n0, n0, n0)
phi(x,y) = 1/(8*pi^2)*cos(2*pi*x)*sin(2*pi*y)
phi0 = function_to_grid(phi, grid)
phi0 = hcat(phi0, phi0, phi0)

# vorticity determined by initial conditions - choose an easy one (parallel laplacian of phi equals vorticity)
w(x,y) = -cos(2*pi*x)*sin(2*pi*y)
w0 = function_to_grid(w, grid)
w0 = hcat(w0, w0, w0)

j(x,y) = 1/(2*pi)*(x^2+y^2-2)*exp(-1/2*(x^2+y^2))
j0 = function_to_grid(j, grid)
j0 = hcat(j0, j0, j0)
psi(x,y) = 1/(2*pi)*exp(-1/2*(x^2+y^2))
psi0 = function_to_grid(psi, grid)
psi0 = hcat(psi0, psi0, psi0)


println(grid.Nk)

Dx = mixed_derivative(1, 0, m, grid)
Dy = mixed_derivative(0, 1, m, grid)

extrall = Extrapolator(M, grid)

Dh = hyperdiffusion(2, mu, grid)
L = Laplacian(m, grid)

wfin, nfin, phifin, psifin, jfin = kelvinhelmholtzcurrent_update!(w0, phi0, n0, psi0, j0, Dx, Dy, extrall, extrall, extrall, extrall, extrall, L, Dh, dt, eta)

anim = @animate for i=1:N
    wfin, nfin, phifin, psifin, jfin = kelvinhelmholtzcurrent_update!(wfin, phifin, nfin, psifin, jfin, Dx, Dy, extrall, extrall, extrall, extrall, extrall, L, Dh, dt, eta)
    heatmap(nfin[:,1], grid)
    println(i)
end

gif(anim, fps=100)