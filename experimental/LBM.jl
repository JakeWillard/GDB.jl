using Plots
using LinearAlgebra
# 1. initialize rho, uvec, f_i and feq_i
# 2. Streaming: move f_i -> f*_i in the direction of evec_i
# 3. Compute macroscopic rho and uvec from f*_i using:
# rho(xvec,t) = sum from i=0:8 of f_i(xvec,t), and
# uvec(xvec,t) = 1/rho * sum from i=0:8 of c*f_i*evec_i
# where c = delta x/delta t is the lattice speed
# 4. compute feq_i using feq_i(xvec,t) = w_i*rho + rho*s_i(uvec(xvec,t))
# 5. calculate updated distribution function f_i = f*_i - 1/tau * (f*_i - feq_i)

# make u, e three dimensional to contain index and 2D vector info?


function lattice_boltzmann_method(xres, yres, f, N, tau, c, count=0)
    # vector at each point with nine dimensions corresponding to each direction 0-8
    
    fnew = zeros((xres, yres, 9))
    for i=1:xres
        for j=1:yres
            # periodic boundary conditions - checking if at corners/edges to ensure no out of bounds
            if i==1 && j==1
                fnew[i,j,1] = f[i,j,1]
                fnew[i,j,2] = f[xres,j,2]
                fnew[i,j,3] = f[i,j+1,3]
                fnew[i,j,4] = f[i+1,j,4]
                fnew[i,j,5] = f[i,yres,5]
                fnew[i,j,6] = f[xres,j+1,6]
                fnew[i,j,7] = f[i+1,j+1,7]
                fnew[i,j,8] = f[i+1,yres,8]
                fnew[i,j,9] = f[xres,yres,9]
            elseif i==1 && j==yres
                fnew[i,j,1] = f[i,j,1]
                fnew[i,j,2] = f[xres,j,2]
                fnew[i,j,3] = f[i,1,3]
                fnew[i,j,4] = f[i+1,j,4]
                fnew[i,j,5] = f[i,j-1,5]
                fnew[i,j,6] = f[xres,1,6]
                fnew[i,j,7] = f[i+1,1,7]
                fnew[i,j,8] = f[i+1,j-1,8]
                fnew[i,j,9] = f[xres,j-1,9]
            elseif i==xres && j==1
                fnew[i,j,1] = f[i,j,1]
                fnew[i,j,2] = f[i-1,j,2]
                fnew[i,j,3] = f[i,j+1,3]
                fnew[i,j,4] = f[1,j,4]
                fnew[i,j,5] = f[i,yres,5]
                fnew[i,j,6] = f[i-1,j+1,6]
                fnew[i,j,7] = f[1,j+1,7]
                fnew[i,j,8] = f[1,yres,8]
                fnew[i,j,9] = f[i-1,yres,9]
            elseif i==xres && j==yres
                fnew[i,j,1] = f[i,j,1]
                fnew[i,j,2] = f[i-1,j,2]
                fnew[i,j,3] = f[i,1,3]
                fnew[i,j,4] = f[1,j,4]
                fnew[i,j,5] = f[i,j-1,5]
                fnew[i,j,6] = f[i-1,1,6]
                fnew[i,j,7] = f[1,1,7]
                fnew[i,j,8] = f[1,j-1,8]
                fnew[i,j,9] = f[i-1,j-1,9]
            elseif i==1
                fnew[i,j,1] = f[i,j,1]
                fnew[i,j,2] = f[xres,j,2]
                fnew[i,j,3] = f[i,j+1,3]
                fnew[i,j,4] = f[i+1,j,4]
                fnew[i,j,5] = f[i,j-1,5]
                fnew[i,j,6] = f[xres,j+1,6]
                fnew[i,j,7] = f[i+1,j+1,7]
                fnew[i,j,8] = f[i+1,j-1,8]
                fnew[i,j,9] = f[xres,j-1,9]
            elseif i==xres
                fnew[i,j,1] = f[i,j,1]
                fnew[i,j,2] = f[i-1,j,2]
                fnew[i,j,3] = f[i,j+1,3]
                fnew[i,j,4] = f[1,j,4]
                fnew[i,j,5] = f[i,j-1,5]
                fnew[i,j,6] = f[i-1,j+1,6]
                fnew[i,j,7] = f[1,j+1,7]
                fnew[i,j,8] = f[1,j-1,8]
                fnew[i,j,9] = f[i-1,j-1,9]
            elseif j==1
                fnew[i,j,1] = f[i,j,1]
                fnew[i,j,2] = f[i-1,j,2]
                fnew[i,j,3] = f[i,j+1,3]
                fnew[i,j,4] = f[i+1,j,4]
                fnew[i,j,5] = f[i,yres,5]
                fnew[i,j,6] = f[i-1,j+1,6]
                fnew[i,j,7] = f[i+1,j+1,7]
                fnew[i,j,8] = f[i+1,yres,8]
                fnew[i,j,9] = f[i-1,yres,9]
            elseif j==yres
                fnew[i,j,1] = f[i,j,1]
                fnew[i,j,2] = f[i-1,j,2]
                fnew[i,j,3] = f[i,1,3]
                fnew[i,j,4] = f[i+1,j,4]
                fnew[i,j,5] = f[i,j-1,5]
                fnew[i,j,6] = f[i-1,1,6]
                fnew[i,j,7] = f[i+1,1,7]
                fnew[i,j,8] = f[i+1,j-1,8]
                fnew[i,j,9] = f[i-1,j-1,9]
            else
                fnew[i,j,1] = f[i,j,1]
                fnew[i,j,2] = f[i-1,j,2]
                fnew[i,j,3] = f[i,j+1,3]
                fnew[i,j,4] = f[i+1,j,4]
                fnew[i,j,5] = f[i,j-1,5]
                fnew[i,j,6] = f[i-1,j+1,6]
                fnew[i,j,7] = f[i+1,j+1,7]
                fnew[i,j,8] = f[i+1,j-1,8]
                fnew[i,j,9] = f[i-1,j-1,9]
            end
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
                u[i,j,k] = c*fnew[i,j,k]/rho[i,j]
            end
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

    s = zeros((xres, yres, 9))
    usquared = zeros((xres, yres))
    for i=1:xres
        for j=1:yres
            for k=1:9
                usquared[i,j] = usquared[i,j] + u[i,j,k]^2
            end
            for m=1:9
                s[i,j,m] = w[m] * (3*u[i,j,m]/c + 9/2*u[i,j,m]^2/c^2 - 3/2*usquared[i,j]/c^2)
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
                fnext[i,j,k] = 1/(1 - 1/tau) * (fnew[i,j,k] - 1/tau*feq[i,j,k])
                # fnext[i,j,k] = fnew[i,j,k] - (1/tau)*(fnew[i,j,k]-feq[i,j,k])
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

xres = Int32(100)
yres = Int32(100)
# gaussian distribution in x and y?
fgauss = zeros(xres, yres, 9)
for i=1:xres
    for j=1:yres
        for k=1:9
            fgauss[i,j,k] = 1/(50^2*sqrt(2*pi)) * exp(-(i-xres/2)^2/(2*50^2)) * exp(-(j-yres/2)^2/(2*50^2))
        end
    end
end
N = 10
tau = .25
deltax = 10000
deltat = 1
c = deltax/deltat

f = lattice_boltzmann_method(xres, yres, fgauss, N, tau, c)
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


anim = @animate for i=1:N
    fgauss = lattice_boltzmann_method(xres, yres, fgauss, 1, tau, c)
    fplot = zeros(xres, yres)
    for i=1:xres
        for j=1:yres
            for k=1:9
                fplot[i,j] = fplot[i,j] + fgauss[i,j,k]
            end
            f0[i,j] = f[i,j,5]
        end
    end
    println(fplot[Int32(xres/2),Int32(yres/2)])
    heatmap(fplot, clims=(0,170))
end
gif(anim, fps=1)