# lbm that calculates eigenvectors/eigenvalues
# cuts off eigenvalues after a certain drop in orders of magnitude (input?)
using Plots
using LinearAlgebra
using SparseArrays
using Arpack

include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/experimental/LBM.jl")

function lattice_boltzmann_basis(xres, yres, f, N, tau, deltax, deltat, threshmag, count=0, rhofull=zeros(xres*yres, N))
    # returns the eigenvalues/eigenvectors that define the basis of the density function over the given time period
    c = deltax/deltat

    ftransform = zeros(xres*yres, 9)
    for j=1:yres
        for i=1:xres
            ftransform[(j-1)*xres+i,:] = f[i,j,:]
        end
    end
    
    stensor = PeriodicStreamTensor(xres, yres)
    fnew = stensor * ftransform

    rho = zeros((xres*yres))
    for i=1:xres*yres
        rho[i] = sum(fnew[i,:])
    end

    u = zeros(xres*yres, 2)
    for i=1:xres*yres
        ux = (fnew[i,2] - fnew[i,4] + 1/sqrt(2) * (fnew[i,6] + fnew[i,9] - fnew[i,7] - fnew[i,8])) / rho[i]
        uy = (fnew[i,3] - fnew[i,5] + 1/sqrt(2) * (fnew[i,6] + fnew[i,8] - fnew[i,9] - fnew[i,8])) / rho[i]
        u[i,1] = ux
        u[i,2] = uy
    end

    unitvecs = zeros(9,2)
    unitvecs[:,1] = [0, c, 0, -c, 0, 1/sqrt(2)*c, -1/sqrt(2)*c, -1/sqrt(2)*c, 1/sqrt(2)*c]
    unitvecs[:,2] = [0, 0, c, 0, -c, 1/sqrt(2)*c, 1/sqrt(2)*c, -1/sqrt(2)*c, -1/sqrt(2)*c]

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
        fnext[i,:] = fnew[i,:] .- (deltat/(tau+.5*deltat)).*(fnew[i,:] - feq[i,:])
    end

    fn = zeros(xres, yres, 9)
    for j=1:yres
        for i=1:xres
            fn[i,j,:] = fnext[(j-1)*xres+i,:]
        end
    end

    if count == 0
        rhofull = rho
    else
        rhofull = hcat(rhofull, rho)
    end

    count = count + 1
    

    if count < N
        lattice_boltzmann_basis(xres, yres, fn, N, tau, deltax, deltat, threshmag, count, rhofull)
    else
        p = size(rhofull)[2]
        covarmatrix = 1/(p-1) * rhofull * transpose(rhofull)
        nevv = 100
        evals, evecs = eigs(covarmatrix, nev=nevv, which=:LM)

        for i=1:nevv-1
            if real(evals[i])/real(evals[1]) < threshmag
                evalsnew = zeros(i)
                evecsnew = zeros(xres*yres,i)
                evalsnew[:] .= evals[1:i,1]
                evecsnew[:,:] .= evecs[:,1:i]
                coeff = transpose(evecsnew) * rhofull
                return evalsnew, evecsnew, coeff
            end
        end
        
    end
end


xres = 100
yres = 100
N = 100
tau = 10000
deltax = 1
deltat = 1
threshmag = 10^-8

fsine = zeros(xres, yres, 9)
for i=1:xres
    for j=1:yres
        for k=1:9
            fsine[i,j,k] = sin(i*pi/xres)*sin(j*pi/yres)
        end
    end
end

evals, evecs, coefficients = lattice_boltzmann_basis(xres, yres, fsine, N, tau, deltax, deltat, threshmag)
println(size(evals))