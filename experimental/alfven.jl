function spline2d(xr, yr, mx, my)

    out = zeros(Float64, mx*my)

    for i=1:mx
        for j=1:my
            k = i + (j-1)*mx
            a = xr^(i-1) / factorial(i-1)
            b = yr^(j-1) / factorial(j-1)
            out[k] = a * b
        end
    end

    return out
end


function stencil2d(mx, my)

    # (i,j) indices for center point
    ic = Int(ceil(mx/2.0))
    jc = Int(ceil(my/2.0))

    M = zeros(Float64, (mx*my, mx*my))

    for i=1:mx
        for j=1:my
            k = i + (j-1)*mx
            M[k,:] = spline2d(i-ic, j-jc, mx, my)
        end
    end

    return Matrix(transpose(inv(M)))
end


function interpolation_row(x, y, mx, my, MinvT, grd::Grid)

    x0, y0 = grd.r0
    dx = grd.dr
    dy = grd.dr
    Nx = grd._Nx
    Ny = grd._Ny

    # (i,j) indices for center point of stencil
    ics = Int(ceil(mx/2.0))
    jcs = Int(ceil(my/2.0))

    # indices for center point in grid
    # NOTE: the addition of 1e-5 before calling floor is to compensate for annoying floating point imprecision.
    icg = Int(floor((x - x0)/dx + 1e-5)) + 1 + grd._Nbuffer
    jcg = Int(floor((y - y0)/dy + 1e-5)) + 1 + grd._Nbuffer

    # values for row
    xr = (x - x0)/dx - (icg - 1 - grd._Nbuffer)
    yr = (y - y0)/dy - (jcg - 1 - grd._Nbuffer)
    spline = spline2d(xr, yr, mx, my)
    row_dat = MinvT * spline

    # collective index for first nonzero element of output
    k0 = (icg - ics + 1) + (jcg - jcs)*Nx

    # indices for nonzero elements of row
    row_js = zeros(Int, mx*my)
    for i=1:mx
        for j=1:my
            row_js[i+(j-1)*mx] = k0 + (i-1) + (j-1)*Nx
        end
    end

    return row_dat, row_js
end



# function bilinear_coefficients(x, y, grd::Grid)
#
#     i = Int64(floor((x - grd.r0[1]) / grd.dr + 1e-5)) + 1 + grd._Nbuffer
#     j = Int64(floor((y - grd.r0[2]) / grd.dr + 1e-5)) + 1 + grd._Nbuffer
#     u = (rem(x - grd.r0[1], grd.dr) + 1e-5) / grd.dr
#     v = (rem(y - grd.r0[2], grd.dr) + 1e-5) / grd.dr
#
#     k0 = i + (j - 1)*grd._Nx
#     ks = Int64[k0, k0+grd._Nx, k0+1, k0+grd._Nx+1]
#     C = [(1 - u)*(1 - v), (1 - u)*v, u*(1 - v), u*v]
#
#     return ks, C
# end


function parallel_laplacian(q::Float64, Nz::Int64, grd::Grid)

    is = Int64[]
    js = Int64[]
    dats = Float64[]

    W = exp(im*2*pi/(q*Nz))
    m = 2
    MinvT = stencil2d(m, m)

    for k=1:grd.Nk

        z = Complex(grd.points[:,k]...)
        z_fwd = W*z
        z_bck = z / W
        ds = abs(z)*2*pi/(q*Nz)

        row_dat, row_js = interpolation_row(real(z), imag(z), m, m, MinvT, grd)
        is = [is; k*ones(Int64, m^2)]
        js = [js; row_js]
        dats = [dats; -2*row_dat / ds^2]

        row_dat, row_js = interpolation_row(real(z_fwd), imag(z_fwd), m, m, MinvT, grd)
        is = [is; k*ones(Int64, m^2)]
        js = [js; row_js]
        dats = [dats; row_dat / ds^2]

        row_dat, row_js = interpolation_row(real(z_bck), imag(z_bck), m, m, MinvT, grd)
        is = [is; k*ones(Int64, m^2)]
        js = [js; row_js]
        dats = [dats; row_dat / ds^2]

        # ks, C = bilinear_coefficients(real(z_fwd), imag(z_fwd), grd)
        # is = [is; k*ones(Int64, 4)]
        # js = [js; ks]
        # dats = [dats; C / ds^2]
        #
        # ks, C = bilinear_coefficients(real(z_bck), imag(z_bck), grd)
        # is = [is; k*ones(Int64, 4)]
        # js = [js; ks]
        # dats = [dats; C / ds^2]

    end

    return sparse(is, js, dats, grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
end


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

        I + (-1)^n*mu*(kron(Iy, Dx) + kron(Dy, Ix))
    end

end


function direct_solve(A, b, bval, extr)

    Anew, bnew = extr(A, b, bval)
    bnew = transpose(A)*b
    Anew = transpose(A)*A
    return extr(Anew \ bnew, bval)
end


function alfven_wave_update!(u, v, Dss, va, Dh, extr, p1, p2)

    bval = zeros(size(u))
    v_t = va^2 * Dss * u
    v[:] = (v + p1.*v_t + p2.*extr(v, bval)) ./ (1 .+ p2)
    v[:] = direct_solve(Dh, v, bval, extr)
    u[:] = (u + p1.*v + p2.*extr(u, bval)) ./ (1 .+ p2)
    u[:] = direct_solve(Dh, u, bval, extr)
end
