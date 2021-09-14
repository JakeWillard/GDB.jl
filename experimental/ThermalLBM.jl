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

function thermal_lattice_boltzmann_method(streamtensor, Dx, Dy, Dxx, Dyy, xres, yres, f, g, tau, tauc, deltax, deltat)
    # vector at each point with nine dimensions corresponding to each direction 0-8
    c = deltax/deltat
    nu = tau*c^2/3
    # xres = grid._Nx
    # yres = grid._Ny

    ftransform = zeros(xres*yres, 9)
    gtransform = zeros(xres*yres, 9)
    for j=1:yres
        for i=1:xres
            ftransform[(j-1)*xres+i,:] = f[i,j,:]
            gtransform[(j-1)*xres+i,:] = g[i,j,:]
        end
    end
    
    
    fnew = streamtensor * ftransform
    gnew = streamtensor * gtransform
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
    # divu = zeros(xres*yres)
    # grad = [0, Dx, Dy, -Dx, -Dy, 1/sqrt(2)*(Dx+Dy), 1/sqrt(2)*(-Dx+Dy), 1/sqrt(2)*(-Dx-Dy), 1/sqrt(2)*(Dx-Dy)]
    # divu = Dx*u[:,2] + Dy*u[:,3] - Dx*u[:,4] - Dy*u[:,5] + 1/sqrt(2).*(Dx*u[:,6]+Dy*u[:,6]) + 1/sqrt(2).*(-Dx*u[:,7]+Dy*u[:,7]) + 1/sqrt(2).*(-Dx*u[:,8]-Dy*u[:,8]) + 1/sqrt(2).*(Dx*u[:,9]-Dy*u[:,9])

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
    """for k=1:9
        esubk = zeros(xres*yres, 9)
        esubk[:,k] .= 1
        for n=1:9
            q[:,k] = q[:,k] + 9*(-(grad[n]*p[:]) ./ rho[:] + nu .* (Dxx*u[:,n]+Dyy*u[:,n])) .* (esubk[:,n]-u[:,n]) + nu .* divu[:] .* (grad[n] * (esubk[:,n]-u[:,n]))
            for m=1:9
                q[:,k] = q[:,k] + (esubk[:,m]-u[:,m]).* (grad[m]*u[:,n]) .*(esubk[:,n]-u[:,n])
            end
        end
    end"""


    
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
        for k=1:9
            usquared[i] = dot(u[i,:],u[i,:])
        end
        for k=1:9
            s[i,k] = w[k] * (3*dot(unitvecs[k,:],u[i,:])/c^2 + 9/2*(dot(unitvecs[k,:],u[i,:]))^2/c^3 - 3/2*usquared[i]/c^2)
        end
    end

    feq = zeros((xres*yres, 9))
    for i=1:xres*yres
        for k=1:9
            feq[i,k] = w[k]*rho[i] + rho[i]*s[i,k]
        end
    end

    fnext = zeros((xres*yres, 9))
    for i=1:xres*yres
        for k=1:9
            # fnext[i,j,k] = 1/(1 - 1/tau) * (fnew[i,j,k] - 1/tau*feq[i,j,k])
            # fnext[i,j,k] = fnew[i,j,k] - (1/tau)*(fnew[i,j,k]-feq[i,j,k])
            fnext[i,k] = fnew[i,k] - (.5*deltat/(tau+.5*deltat))*(fnew[i,k]-feq[i,k])
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
        for k=1:9
            gnext[i,k] = gnew[i,k] - deltat/(tauc + .5*deltat)*(gnew[i,k] - geq[i,k]) - tauc/(tauc + .5*deltat)*fnext[i,k]*q[i,k]*deltat
        end
    end

    fn = zeros(xres, yres, 9)
    gn = zeros(xres, yres, 9)
    ex = zeros(xres, yres, 9)
    ex1 = zeros(xres, yres)
    ex2 = zeros(xres, yres, 2)
    for j=1:yres
        for i=1:xres
            fn[i,j,:] = fnext[(j-1)*xres+i,:]
            gn[i,j,:] = gnext[(j-1)*xres+i,:]
            ex[i,j,:] = q[(j-1)*xres+i,:]
            ex1[i,j] = rho[(j-1)*xres+i]
            ex2[i,j,:] = u[(j-1)*xres+i,:]
        end
    end
    

    return fn, gn, ex, ex1, ex2
end


N = 200
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

stensor = MirrorStreamTensor(xres, yres)

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
# heatmap(gsine)

anim = @animate for i=1:N
    fsine, gsine, ex, ex1, ex2 = thermal_lattice_boltzmann_method(stensor, pDx, pDy, pDxx, pDyy, xres, yres, fsine, gsine, tau, tauc, deltax, deltat)
    fplot = zeros(xres, yres)
    gplot = zeros(xres, yres)
    explot = zeros(xres, yres)
    exp1 = zeros(xres, yres)
    exp2 = zeros(xres, yres)
    ucheck = zeros(xres, yres)
    for i=1:xres
        for j=1:yres
            exp1[i,j] = ex1[i,j]
            # explot[i,j] = ex[i,j,1]
            exp2[i,j] = ex2[i,j,2]
            for k=1:9
                fplot[i,j] = fplot[i,j] + fsine[i,j,k]
                gplot[i,j] = gplot[i,j] + gsine[i,j,k]
                explot[i,j] = explot[i,j] + ex[i,j,k]
                ucheck[i,j] = ucheck[i,j] + ex[i,j,k]/ex1[i,j]
            end
        end
    end
    println(fplot[Int32(xres/2),Int32(yres/2)])
    println(gplot[Int32(xres/2),Int32(yres/2)])
    heatmap(gplot, clims=(0,8))
end
gif(anim, fps=10)