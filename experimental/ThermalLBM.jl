using SparseArrays
using LinearAlgebra
using Plots
using RecipesBase


include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/domain/grid.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/domain/barrier.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/domain/operators.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/domain/interpolation.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/experimental/LBM.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/experimental/periodic_derivatives.jl")
include("C:/Users/lucas/OneDrive/Documents/Dartmouth/Research/Code/LBMEigen.jl")

function thermal_lattice_boltzmann_method(streamtensor, Dx, Dy, Dxx, Dyy, xres, yres, f, g, tau, tauc, deltax, deltat)
    # vector at each point with nine dimensions corresponding to each direction 0-8
    c = deltax/deltat
    nu = tau*c^2/3
    # xres = grid._Nx
    # yres = grid._Ny
    
    
    fnew = streamtensor * f
    gnew = streamtensor * g
    # println(size(ftransform))


    rho = zeros((xres*yres))
    for i=1:xres*yres
        rho[i] = sum(fnew[i,:])
    end

    p = zeros(xres*yres)
    for i=1:xres*yres
        p[i] = c^2/3 * rho[i]
    end

    u = zeros(xres*yres, 2)
    for i=1:xres*yres
        ux = (fnew[i,2] - fnew[i,4] + 1/sqrt(2) * (fnew[i,6] + fnew[i,9] - fnew[i,7] - fnew[i,8])) / rho[i]
        uy = (fnew[i,3] - fnew[i,5] + 1/sqrt(2) * (fnew[i,6] + fnew[i,8] - fnew[i,9] - fnew[i,8])) / rho[i]
        u[i,1] = ux
        u[i,2] = uy
    end
    
    # derivative matrices are operators that act on other matrices
    
    q = zeros(xres*yres, 9)

    grad = [Dx, Dy]
    divu = Dx*u[:,1] + Dy*u[:,2]
    gradp = zeros(xres*yres, 2)
    gradp[:,1] = Dx*p[:] 
    gradp[:,2] = Dy*p[:]
    laplu = zeros(xres*yres, 2)
    laplu[:,1] = (Dxx+Dyy)*u[:,1]
    laplu[:,2] = (Dxx+Dyy)*u[:,2]

    # unitvecs = [[0,0] [c,0] [0,c] [-c,0] [0,-c] 1/sqrt(2)*[c,c] 1/sqrt(2)*[-c,c] 1/sqrt(2)*[-c,-c] 1/sqrt(2)*[c,-c]]
    unitvecs = zeros(9,2)
    unitvecs[:,1] = [0, c, 0, -c, 0, 1/sqrt(2)*c, -1/sqrt(2)*c, -1/sqrt(2)*c, 1/sqrt(2)*c]
    unitvecs[:,2] = [0, 0, c, 0, -c, 1/sqrt(2)*c, 1/sqrt(2)*c, -1/sqrt(2)*c, -1/sqrt(2)*c]

    vgradu = zeros(xres*yres, 9, 2)
    for k=1:9
        vgradu[:,k,1] = (unitvecs[k,1] .- u[:,1]) .* (Dx*u[:,1]) + (unitvecs[k,2] .- u[:,2]) .* (Dy*u[:,1])
        vgradu[:,k,2] = (unitvecs[k,1] .- u[:,1]) .* (Dx*u[:,2]) + (unitvecs[k,2] .- u[:,2]) .* (Dy*u[:,2])
    end
    

    # .* for element wise multiplication, similar for division
    for i=1:xres*yres
        for a=1:9
            q[i,a] = - 1/rho[i]*dot(unitvecs[a,:]-u[i,:], gradp[i,:]) 
            + nu*(dot(unitvecs[a,:]-u[i,:], laplu[i,:]))
            + dot(unitvecs[a,:]-u[i,:], vgradu[i,a,:])
        end
    end

    w = zeros(9)
    w[1] = 4/9
    w[2] = 1/9
    w[3] = 1/9
    w[4] = 1/9
    w[5] = 1/9
    w[6] = 1/36
    w[7] = 1/36
    w[8] = 1/36
    w[9] = 1/36

    s = zeros((xres*yres, 9))
    usquared = zeros((xres*yres))
    for i=1:xres*yres
        usquared[i] = dot(u[i,:],u[i,:])
        s[i,:] = w[:] .* (3*(unitvecs*u[i,:])/c^2 .+ 9/2*((unitvecs*u[i,:])).^2/c^3 .- 3/2*usquared[i]/c^2)
    end

    feq = zeros((xres*yres, 9))
    for i=1:xres*yres
        feq[i,:] = w[:].*rho[i] .+ rho[i].*s[i,:]
    end

    fnext = zeros((xres*yres, 9))
    for i=1:xres*yres
        for k=1:9
            # fnext[i,j,k] = 1/(1 - 1/tau) * (fnew[i,j,k] - 1/tau*feq[i,j,k])
            # fnext[i,j,k] = fnew[i,j,k] - (1/tau)*(fnew[i,j,k]-feq[i,j,k])
            fnext[i,:] = fnew[i,:] .- (.5*deltat/(tau+.5*deltat)).*(fnew[i,:].-feq[i,:])
        end
    end

    rhoepsilon = zeros(xres*yres)
    for i=1:xres*yres
        rhoepsilon[i] = sum(gnew[i,:] .- deltat/2 .*fnext[i,:].*q[i,:])
    end

    geq = zeros((xres*yres, 9))
    for i=1:xres*yres
        geq[i,1] = -2/3*rhoepsilon[i] * usquared[i]/c^2
        geq[i,2] = rhoepsilon[i]/9 * (1.5 + 1.5*dot(unitvecs[2,:],u[i,:])/c^2 + 4.5*(dot(unitvecs[2,:],u[i,:]))^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,3] = rhoepsilon[i]/9 * (1.5 + 1.5*dot(unitvecs[3,:],u[i,:])/c^2 + 4.5*(dot(unitvecs[3,:],u[i,:]))^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,4] = rhoepsilon[i]/9 * (1.5 + 1.5*dot(unitvecs[4,:],u[i,:])/c^2 + 4.5*(dot(unitvecs[4,:],u[i,:]))^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,5] = rhoepsilon[i]/9 * (1.5 + 1.5*dot(unitvecs[5,:],u[i,:])/c^2 + 4.5*(dot(unitvecs[5,:],u[i,:]))^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,6] = rhoepsilon[i]/36 * (3 + 6*dot(unitvecs[6,:],u[i,:])/c^2 + 4.5*(dot(unitvecs[6,:],u[i,:]))^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,7] = rhoepsilon[i]/36 * (3 + 6*dot(unitvecs[7,:],u[i,:])/c^2 + 4.5*(dot(unitvecs[7,:],u[i,:]))^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,8] = rhoepsilon[i]/36 * (3 + 6*dot(unitvecs[8,:],u[i,:])/c^2 + 4.5*(dot(unitvecs[8,:],u[i,:]))^2/c^4 - 1.5*usquared[i]/c^2)
        geq[i,9] = rhoepsilon[i]/36 * (3 + 6*dot(unitvecs[9,:],u[i,:])/c^2 + 4.5*(dot(unitvecs[9,:],u[i,:]))^2/c^4 - 1.5*usquared[i]/c^2)
    end

    gnext = zeros(xres*yres, 9)
    for i=1:xres*yres
        gnext[i,:] = gnew[i,:] .- deltat/(tauc + .5*deltat).*(gnew[i,:] .- geq[i,:]) .- tauc/(tauc + .5*deltat)*deltat.*fnext[i,:].*q[i,:]
    end
    

    return fnext, gnext, q, rho, usquared
end


N = 100
tau = 1000
tauc = 1000
deltax = 1
deltat = 1

grd = Grid(100, 100, 1; Nbuffer=0) do 

    r0 = Float64[-1,-1]
    r1 = Float64[1,1]
    inside(x,y) = true
    inside, r0, r1
    
end

# xres = grd._Nx
# yres = grd._Ny

xres = 100
yres = 100

stensor = PeriodicStreamTensor(xres, yres)

sx = 5
sy = 5
dx = 1
dy = 1
pDx = periodic_derivative(1, 0, sx, sy, xres, yres, dx, dy)
pDy = periodic_derivative(0, 1, sx, sy, xres, yres, dx, dy)
pDxx = periodic_derivative(2, 0, sx, sy, xres, yres, dx, dy)
pDyy = periodic_derivative(0, 2, sx, sy, xres, yres, dx, dy)

Minvx = stencil1d(sx)
Minvy = stencil1d(sy)
Dx = derivative_matrix(1, 0, Minvx, Minvy, grd)
Dy = derivative_matrix(0, 1, Minvx, Minvy, grd)
Dxx = derivative_matrix(2, 0, Minvx, Minvy, grd)
Dyy = derivative_matrix(0, 2, Minvx, Minvy, grd)

fsine = zeros(xres, yres, 9)
gsine = zeros(xres, yres, 9)
rhoinit = zeros(xres, yres)
uinit = zeros(xres, yres, 9)
for i=1:xres
    for j=1:yres
        for k=1:9
            fsine[i,j,k] = .5 + .4*sin(i*2*pi/(xres))*sin(j*2*pi/(yres))
        end
        rhoinit[i,j] = sum(fsine[i,j,:])
        uinit[i,j,:] = fsine[i,j,:] ./ rhoinit[i,j]
        for k=1:9
            gsine[i,j,k] = (1 - uinit[i,j,k])^2/2 * fsine[i,j,k]
        end
    end
end

fs2 = zeros(xres*yres, 9)
gs2 = zeros(xres*yres, 9)

for j=1:yres
    for i=1:xres
        fs2[(j-1)*xres+i,:] = fsine[i,j,:]
        gs2[(j-1)*xres+i,:] = gsine[i,j,:]
    end
end

rhofull = zeros(xres*yres)
anim = @animate for i=1:N
    fs2, gs2, ex, ex1, ex2 = thermal_lattice_boltzmann_method(stensor, pDx, pDy, pDxx, pDyy, xres, yres, fs2, gs2, tau, tauc, deltax, deltat)
    fplot = zeros(xres, yres)
    gplot = zeros(xres, yres)
    explot = zeros(xres, yres)
    exp1 = zeros(xres, yres)
    exp2 = zeros(xres, yres)
    ucheck = zeros(xres, yres)
    if i==1
        rhofull = ex2
    else
        rhofull = hcat(rhofull, ex2)
    end
    for i=1:xres
        for j=1:yres
            exp1[i,j] = ex1[(j-1)*xres+i]
            exp2[i,j] = ex2[(j-1)*xres+i]
            for k=1:9
                fplot[i,j] = fplot[i,j] + fs2[(j-1)*xres+i,k]
                gplot[i,j] = gplot[i,j] + gs2[(j-1)*xres+i,k]
                explot[i,j] = explot[i,j] + ex[(j-1)*xres+i,k]
                ucheck[i,j] = ucheck[i,j] + ex[(j-1)*xres+i,k]/ex1[(j-1)*xres+i]
            end
        end
    end
    println(fplot[Int32(xres/2),Int32(yres/2)])
    println(gplot[Int32(xres/2),Int32(yres/2)])
    heatmap(gplot, clims=(0,5))
end
# gif(anim, fps=10)

p = size(rhofull)[2]
covarmatrix = 1/(p-1) * rhofull * transpose(rhofull)
evals, evecs = eigs(covarmatrix, nev=20, which=:LM)

evecstransform = zeros(xres, yres)
for j=1:yres
    for i=1:xres
        evecstransform[i,j] = evecs[(j-1)*xres+i,3]
    end
end

heatmap(evecstransform)