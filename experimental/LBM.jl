using Plots
using LinearAlgebra
using SparseArrays
# 1. initialize rho, uvec, f_i and feq_i
# 2. Streaming: move f_i -> f*_i in the direction of evec_i
# 3. Compute macroscopic rho and uvec from f*_i using:
# rho(xvec,t) = sum from i=0:8 of f_i(xvec,t), and
# uvec(xvec,t) = 1/rho * sum from i=0:8 of c*f_i*evec_i
# where c = delta x/delta t is the lattice speed
# 4. compute feq_i using feq_i(xvec,t) = w_i*rho + rho*s_i(uvec(xvec,t))
# 5. calculate updated distribution function f_i = f*_i - 1/tau * (f*_i - feq_i)

# make u, e three dimensional to contain index and 2D vector info?

struct StreamTensor

    E :: SparseMatrixCSC
    SE :: SparseMatrixCSC
    S :: SparseMatrixCSC
    SW :: SparseMatrixCSC
    W :: SparseMatrixCSC
    NW :: SparseMatrixCSC
    N :: SparseMatrixCSC
    NE :: SparseMatrixCSC

end


function PeriodicStreamTensor(Nx, Ny)

    Fx = spdiagm(-1 => ones(Nx - 1))
    Fx[1,Nx] = 1
    Bx = spdiagm(1 => ones(Nx - 1))
    Bx[Nx,1] = 1

    Fy = spdiagm(-1 => ones(Ny - 1))
    Fy[1,Ny] = 1
    By = spdiagm(1 => ones(Ny - 1))
    By[Ny,1] = 1

    Ix = sparse(I, Nx, Nx)
    Iy = sparse(I, Ny, Ny)

    E = kron(Iy, Fx)
    SE = kron(By, Fx)
    S = kron(By, Ix)
    SW = kron(By, Bx)
    W = kron(Iy, Bx)
    NW = kron(Fy, Bx)
    N = kron(Fy, Ix)
    NE = kron(Fy, Fx)

    return StreamTensor(E, SE, S, SW, W, NW, N, NE)
end


function MirrorStreamTensor(Nx, Ny)

    Fx = spdiagm(-1 => ones(Nx - 1))
    Fx[1,1] = 1
    Bx = spdiagm(1 => ones(Nx - 1))
    Bx[Nx,Nx] = 1

    Fy = spdiagm(-1 => ones(Ny - 1))
    Fy[1,1] = 1
    By = spdiagm(1 => ones(Ny - 1))
    By[Ny,Ny] = 1

    Ix = sparse(I, Nx, Nx)
    Iy = sparse(I, Ny, Ny)

    E = kron(Iy, Fx)
    SE = kron(By, Fx)
    S = kron(By, Ix)
    SW = kron(By, Bx)
    W = kron(Iy, Bx)
    NW = kron(Fy, Bx)
    N = kron(Fy, Ix)
    NE = kron(Fy, Fx)

    return StreamTensor(E, SE, S, SW, W, NW, N, NE)
end


function Base.:*(S::StreamTensor, f::Matrix{Float64})

    out = zeros(size(f))

    out[:,1] = f[:,1]
    out[:,2] = S.E * f[:,2]
    out[:,3] = S.N * f[:,3]
    out[:,4] = S.W * f[:,4]
    out[:,5] = S.S * f[:,5]
    out[:,6] = S.NE * f[:,6]
    out[:,7] = S.NW * f[:,7]
    out[:,8] = S.SW * f[:,8]
    out[:,9] = S.SE * f[:,9]

    return out
end


function lattice_boltzmann_method(xres, yres, f, N, tau, deltax, deltat, count=0)
    # vector at each point with nine dimensions corresponding to each direction 0-8
    c = deltax/deltat

    ftransform = zeros(xres*yres, 9)
    for j=1:yres
        for i=1:xres
            ftransform[(j-1)*xres+i,:] = f[i,j,:]
        end
    end
    
    stensor = PeriodicStreamTensor(xres, yres)
    fnewtransform = stensor * ftransform

    fnew = zeros(xres, yres, 9)
    for j=1:yres
        for i=1:xres
            fnew[i,j,:] = fnewtransform[(j-1)*xres+i,:]
        end
    end

    rho = zeros((xres, yres))
    for i=1:xres
        for j=1:yres
            for k=1:9
                rho[i,j] = rho[i,j] + fnew[i,j,k]
            end
        end
    end

    u = zeros((xres, yres, 9))
    for i=1:xres
        for j=1:yres
            for k=1:9
                u[i,j,k] = fnew[i,j,k]/rho[i,j]
                # u[i,j,k] = c*fnew[i,j,k]/rho[i,j]
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

    """w[1] = 4/9
    w[2] = 1/9
    w[3] = 1/9
    w[4] = 1/9
    w[5] = 1/9
    w[6] = 1/36
    w[7] = 1/36
    w[8] = 1/36
    w[9] = 1/36"""

    s = zeros((xres, yres, 9))
    usquared = zeros((xres, yres))
    for i=1:xres
        for j=1:yres
            for k=1:9
                usquared[i,j] = usquared[i,j] + u[i,j,k]^2
            end
            for k=1:9
                s[i,j,k] = w[k] * (3*u[i,j,k]/c + 9/2*u[i,j,k]^2/c^2 - 3/2*usquared[i,j]/c^2)
            end
        end
    end


    feq = zeros((xres, yres, 9))
    for i=1:xres
        for j=1:yres
            for k=1:9
                feq[i,j,k] = w[k]*rho[i,j] + rho[i,j]*s[i,j,k]
            end
        end
    end

    fnext = zeros((xres, yres, 9))
    for i=1:xres
        for j=1:yres
            for k=1:9
                # fnext[i,j,k] = 1/(1 - 1/tau) * (fnew[i,j,k] - 1/tau*feq[i,j,k])
                # fnext[i,j,k] = fnew[i,j,k] - (1/tau)*(fnew[i,j,k]-feq[i,j,k])
                fnext[i,j,k] = fnew[i,j,k] - (deltat/(tau+.5*deltat))*(fnew[i,j,k]-feq[i,j,k])
            end
        end
    end

    count = count + 1
    if count < N
        lattice_boltzmann_method(xres, yres, fnext, N, tau, c, count)
    else
        return fnext
    end
end

xres = Int32(200)
yres = Int32(200)
# gaussian distribution in x and y?
fgauss = zeros(xres, yres, 9)
for i=1:xres
    for j=1:yres
        for k=1:9
            fgauss[i,j,k] = 1/(50^2*sqrt(2*pi)) * exp(-(i-xres/2)^2/(2*50^2)) * exp(-(j-yres/2)^2/(2*50^2))
        end
    end
end

fsine = zeros(xres, yres, 9)
for i=1:xres
    for j=1:yres
        for k=1:9
            fsine[i,j,k] = sin(i*pi/xres)*sin(j*pi/yres)
        end
    end
end

N = 200
tau = 10000
deltax = 1
deltat = 1
c = deltax/deltat

"""f = lattice_boltzmann_method(xres, yres, fgauss, N, tau, deltax, deltat)
fplot = zeros(xres, yres)
f0 = zeros(xres, yres)
for i=1:xres
    for j=1:yres
        for k=1:9
            fplot[i,j] = fplot[i,j] + f[i,j,k]
        end
        f0[i,j] = f[i,j,5]
    end
end
"""

"""anim = @animate for i=1:N
    fsine = lattice_boltzmann_method(xres, yres, fsine, 1, tau, deltax, deltat)
    fplot = zeros(xres, yres)
    for i=1:xres
        for j=1:yres
            for k=1:9
                fplot[i,j] = fplot[i,j] + fsine[i,j,k]
            end
            f0[i,j] = f[i,j,5]
        end
    end
    println(fplot[Int32(xres/2),Int32(yres/2)])
    heatmap(fplot)
end
gif(anim, fps=10)"""