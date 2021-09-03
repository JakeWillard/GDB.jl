using SparseArrays
using LinearAlgebra
using Plots
using RecipesBase


include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/domain/grid.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/domain/barrier.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/domain/operators.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/domain/interpolation.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/experimental/LBM.jl")


function thermal_lattice_boltzmann_method(grid::Grid, f, g, p, N, tau, tauc, deltax, deltat, count=0)
    # vector at each point with nine dimensions corresponding to each direction 0-8
    c = deltax/deltat
    nu = tau*c^2/3
    xres = grid._Nx
    yres = grid._Ny

    ftransform = zeros(xres*yres, 9)
    gtransform = zeros(xres*yres, 9)
    for j=1:yres
        for i=1:xres
            ftransform[(j-1)*xres+i,:] = f[i,j,:]
            gtransform[(j-1)*xres+i,:] = g[i,j,:]
        end
    end
    
    stensor = PeriodicStreamTensor(xres, yres)
    fnew = stensor * ftransform
    gnew = stensor * gtransform
    # println(size(ftransform))


    rho = zeros((xres*yres))
    for i=1:xres*yres
        rho[i] = sum(fnew[i,:])
    end

    u = zeros((xres*yres, 9))
    for i=1:xres
        for k=1:9
            u[i,k] = fnew[i,k]./rho[i]
            # u[i,j,k] = c*fnew[i,j,k]/rho[i,j]
        end
    end
    
    Minvx = stencil1d(5)
    Minvy = stencil1d(5)
    spDx = derivative_matrix(1, 0, Minvx, Minvy, grid)
    spDy = derivative_matrix(0, 1, Minvx, Minvy, grid)
    spDxx = derivative_matrix(2, 0, Minvx, Minvy, grid)
    spDyy = derivative_matrix(0, 2, Minvx, Minvy, grid)
    Dx = Matrix(spDx)
    Dy = Matrix(spDy)
    Dxx = Matrix(spDxx)
    Dyy = Matrix(spDyy)
    # derivative matrices are operators that act on other matrices

    ptransform = zeros(xres*yres)
    for j=1:yres
        for i=1:xres
            ptransform[(j-1)*xres+i] = p[i,j]
        end
    end
    
    q = zeros(xres*yres, 9)
    # divu = zeros(xres*yres)
    grad = [0, Dx, Dy, -Dx, -Dy, 1/sqrt(2)*(Dx+Dy), 1/sqrt(2)*(-Dx+Dy), 1/sqrt(2)*(-Dx-Dy), 1/sqrt(2)*(Dx-Dy)]
    divu = Dx*u[:,2] + Dy*u[:,3] - Dx*u[:,4] - Dy*u[:,5] + 1/sqrt(2).*(Dx*u[:,6]+Dy*u[:,6]) + 1/sqrt(2).*(-Dx*u[:,7]+Dy*u[:,7]) + 1/sqrt(2).*(-Dx*u[:,8]-Dy*u[:,8]) + 1/sqrt(2).*(Dx*u[:,9]-Dy*u[:,9])

    # .* for element wise multiplication, similar for division
    for k=1:9
        esubk = zeros(xres*yres, 9)
        esubk[:,k] .= 1
        for n=1:9
            q[:,k] = q[:,k] + 9*((grad[n]*ptransform[:]) ./ rho[:] + nu .* (Dxx*u[:,n]+Dyy*u[:,n])) .* (esubk[:,n]-u[:,n]) + nu .* divu[:] .* (grad[n] * (esubk[:,n]-u[:,n]))
            for m=1:9
                q[:,k] = q[:,k] + (esubk[:,m]-u[:,m]).* (grad[m]*u[:,n]) .*(esubk[:,n]-u[:,n])
            end
        end
    end

    w = zeros(9)
    w[1] = 0
    w[2] = c
    w[3] = c
    w[4] = c
    w[5] = c
    w[6] = c
    w[7] = c
    w[8] = c
    w[9] = c

    s = zeros((xres*yres, 9))
    usquared = zeros((xres*yres))
    for i=1:xres*yres
        for k=1:9
            usquared[i] = usquared[i] + u[i,k]^2
        end
        for k=1:9
            s[i,k] = w[k] * (3*u[i,k]/c^2 + 9/2*u[i,k]^2/c^4 - 3/2*usquared[i]/c^2)
        end
    end

    feq = zeros((xres*yres, 9))
    for i=1:xres*yres
        for k=1:9
            feq[i,k] = w[k]*rho[i] + rho[i].*s[i,k]
        end
    end

    fnext = zeros((xres*yres, 9))
    for i=1:xres*yres
        for k=1:9
            # fnext[i,j,k] = 1/(1 - 1/tau) * (fnew[i,j,k] - 1/tau*feq[i,j,k])
            # fnext[i,j,k] = fnew[i,j,k] - (1/tau)*(fnew[i,j,k]-feq[i,j,k])
            fnext[i,k] = fnew[i,k] - (deltat/(tau+.5*deltat))*(fnew[i,k]-feq[i,k])
        end
    end

    rhoepsilon = zeros(xres*yres)
    for i=1:xres*yres
        for k=1:9
            rhoepsilon[i] = rhoepsilon[i] + (gnew[i,k] - deltat/2*fnext[i,k]*q[i,k])
        end
    end

    geq = zeros((xres*yres, 9))
    for i=1:xres*yres
        geq[i,1] = -2/3*rhoepsilon[i] * usquared[i]/c^2
        geq[i,2] = rhoepsilon[i]/9 * (1.5 + 1.5*u[i,2]/c^2 + 4.5*(u[i,2])^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,3] = rhoepsilon[i]/9 * (1.5 + 1.5*u[i,3]/c^2 + 4.5*(u[i,3])^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,4] = rhoepsilon[i]/9 * (1.5 + 1.5*u[i,4]/c^2 + 4.5*(u[i,4])^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,5] = rhoepsilon[i]/9 * (1.5 + 1.5*u[i,5]/c^2 + 4.5*(u[i,5])^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,6] = rhoepsilon[i]/36 * (3 + 6*u[i,6]/c^2 + 4.5*(u[i,6])^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,7] = rhoepsilon[i]/36 * (3 + 6*u[i,7]/c^2 + 4.5*(u[i,7])^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,8] = rhoepsilon[i]/36 * (3 + 6*u[i,8]/c^2 + 4.5*(u[i,8])^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,9] = rhoepsilon[i]/36 * (3 + 6*u[i,9]/c^2 + 4.5*(u[i,9])^2/c^4 - 1.5*usquared[i]/c^2)
    end

    gnext = zeros(xres*yres, 9)
    for i=1:xres*yres
        for k=1:9
            gnext[i,k] = gnew[i,k] - deltat/(tauc + .5*deltat)*(gnew[i,k] - geq[i,k]) - tauc/(tauc + .5*deltat)*fnext[i,k]*q[i,k]*deltat
        end
    end

    fn = zeros(xres, yres, 9)
    gn = zeros(xres, yres, 9)
    for j=1:yres
        for i=1:xres
            fn[i,j,:] = fnext[(j-1)*xres+i,:]
            gn[i,j,:] = gnext[(j-1)*xres+i,:]
        end
    end
    

    count = count + 1
    if count < N
        thermal_lattice_boltzmann_method(grid, fn, gn, p, N, tau, tauc, deltax, deltat, count)
    else
        return fn, gn
    end
end

N = 10
tau = 1000
tauc = 1000
deltax = 1000
deltat = 1

grd = Grid(100, 100, 1; Nbuffer=0) do 

    r0 = Float64[-1,-1]
    r1 = Float64[1,1]
    inside(x,y) = true
    inside, r0, r1
    
end

xres = grd._Nx
yres = grd._Ny

fsine = zeros(xres, yres, 9)
gsine = zeros(xres, yres, 9)
p = zeros(xres, yres)
for i=1:xres
    for j=1:yres
        p[i,j] = sin(i*pi/xres)*sin(j*pi/yres)
        for k=1:9
            fsine[i,j,k] = sin(i*pi/xres)*sin(j*pi/yres)
            gsine[i,j,k] = sin(i*pi/xres)*sin(j*pi/yres)
        end
    end
end

anim = @animate for i=1:N
    fsine, gsine = thermal_lattice_boltzmann_method(grd, fsine, gsine, p, 1, tau, tauc, deltax, deltat)
    fplot = zeros(xres, yres)
    gplot = zeros(xres, yres)
    for i=1:xres
        for j=1:yres
            for k=1:9
                fplot[i,j] = fplot[i,j] + fsine[i,j,k]
                gplot[i,j] = gplot[i,j] + gsine[i,j,k]
            end
        end
    end
    println(fplot[Int32(xres/2),Int32(yres/2)])
    println(gplot[Int32(xres/2),Int32(yres/2)])
    heatmap(gplot)
end
gif(anim, fps=10)