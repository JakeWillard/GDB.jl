


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


function hyperdiffusion(n, mu, grd)

    operator_to_grid(grd) do
        Dx = derivative_1d(grd._Nx, 2*n, 2*n+1)
        Dy = derivative_1d(grd._Ny, 2*n, 2*n+1)
        Ix = sparse(I, grd._Nx, grd._Nx)
        Iy = sparse(I, grd._Ny, grd._Ny)

        I - mu*(kron(Iy, Dx) + kron(Dy, Ix))
    end

end


function direct_solve(A, b, bval, extr)

    Anew, bnew = extr(A, b, bval)
    bnew = transpose(A)*b
    Anew = transpose(A)*A
    return extr(Anew \ bnew, bval)
end


# function wave_update!(u, Dxx, Dyy, cx, cy, Dh, extr)
#
#     bval = zeros(size(u)[1])
#     u[:,3] = 2*u[:,2] - u[:,1] + (cx^2*Dxx + cy^2*Dyy)*u[:,2]
#     # u[:,3] = direct_solve(Dh, u[:,3], bval, extr)
#     u[:,1:2] = u[:,2:3]
# end


function wave_update!(u, v, Dxx, Dyy, cx, cy, Dh, extr)

    bval = zeros(size(u))
    v_t = (cx^2*Dxx + cy^2*Dyy) * u
    v[:] = (v + p1.*v_t + p2.*extr(v, bval)) ./ (1 .+ p2)
    v[:] = direct_solve(Dh, v, bval, extr)
    u[:] = (u + p1.*v + p2.*extr(u, bval)) ./ (1 .+ p2)
    u[:] = direct_solve(Dh, u, bval, extr)
end



function wakatani_update!(w, phi, n, Dx, Dy, L, Dh, w_ex, n_ex, phi_ex, alpha, kx, ky, dt)

    bval = zeros(size(w)[1])

    phi_x = Dx*phi[:,2]
    phi_y = Dy*phi[:,2]
    n_x = Dx*n[:,2]
    n_y = Dy*n[:,2]
    w_x = Dx*w[:,2]
    w_y = Dy*w[:,2]

    w_t = -(phi_x.*w_y - phi_y.*w_x) + alpha*(phi - n[:,2])
    n_t = -(phi_x.*n_y - phi_y.*n_x) - (ky*phi_x - kx*phi_y) + alpha*(phi - n[:,2])

    w[:,3] = (w[:,1] + w[:,2])/2 + dt*w_t
    n[:,3] = (n[:,1] + n[:,2])/2 + dt*n_t

    w[:,3] = direct_solve(Dh, w[:,3], bval, w_ex)
    n[:,3] = direct_solve(Dh, n[:,3], bval, n_ex)
    phi[:] = direct_solve(L, w[:,3], bval, phi_ex)

    phi_x = Dx*phi[:,3]
    phi_y = Dy*phi[:,3]
    n_x = Dx*n[:,3]
    n_y = Dy*n[:,3]
    w_x = Dx*w[:,3]
    w_y = Dy*w[:,3]

    w_t = -(phi_x.*w_y - phi_y.*w_x) + alpha*(phi - n[:,3])
    n_t = -(phi_x.*n_y - phi_y.*n_x) - (ky*phi_x - kx*phi_y) + alpha*(phi - n[:,3])

    w[:,3] = w[:,2] + dt*w_t
    n[:,3] = n[:,2] + dt*n_t

    w[:,3] = direct_solve(Dh, w[:,3], bval, w_ex)
    n[:,3] = direct_solve(Dh, n[:,3], bval, n_ex)
    phi[:] = direct_solve(L, w[:,3], bval, phi_ex)

    w[:,1:2] = w[:,2:3]
    n[:,1:2] = n[:,2:3]
end


function linsolve_timing(N, m, alpha, fraction, Np; Nsolves=10)

    A = spdiagm(0 => rand(N))
    for i=1:m
        A += spdiagm(i => rand(N - i)) / alpha^i
        A += spdiagm(-i => rand(N - i)) / alpha^i
    end
    println("created matrix")

    x0 = zeros(N)
    b = rand(N)

    stat = @timed begin
        for _=1:Nsolves
            x = prbgs(A, x0, b; block_fraction=fraction, Np=Np)
        end
    end

    return stat.time
end


function big_matrix_with_pmap(N, m, alpha, Nit)

    A = spdiagm(0 => rand(N))
    for i=1:m
        A += spdiagm(i => rand(N - i)) / alpha^i
        A += spdiagm(-i => rand(N - i)) / alpha^i
    end
    B = kron(A, A)

    x = pmap(1:Nit) do i
        2*B
        i
    end
end


function big_matrix_with_darray(N, m, alpha, Nit)

    A = spdiagm(0 => rand(N))
    for i=1:m
        A += spdiagm(i => rand(N - i)) / alpha^i
        A += spdiagm(-i => rand(N - i)) / alpha^i
    end
    B = kron(A, A)

    x = DArray((Nit,), workers(), length(workers())) do inds
        2*A
        [inds...]
    end
end
