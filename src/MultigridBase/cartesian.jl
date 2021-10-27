
function cartesian_1d_crs(N, m)

    # make stencil
    ic = Int(ceil(m/2.0))
    M = zeros((m, m))
    for i=1:m
        for j=1:m
            M[i, j] = (i - ic)^(j-1) / factorial(j-1)
        end
    end
    Minv = inv(M)
    sten = transpose(Minv) * Float64[(0.5)^(j-1)/factorial(j-1) for j=1:m]

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
    Rout = transpose(Iout) / 2

    return sparse(Iout), sparse(Rout)
end


function cartesian_2d_crs(Nx, Ny, mx, my)

    Ix, Rx = cartesian_1d_crs(Nx, mx)
    Iy, Ry = cartesian_1d_crs(Ny, my)

    I2d = kron(Iy, Ix)
    R2d = kron(Ry, Rx)
    return I2d, R2d
end
