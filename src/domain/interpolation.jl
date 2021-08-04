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
    dx = grd.dx
    dy = grd.dy
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
    xr = (x - x0)/dx - (icg - 1)
    yr = (y - y0)/dy - (jcg - 1)
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


function interpolation_function(V::Vector{Float64}, mx, my, MinvT, grd::Grid)

    return (x, y) -> dot(V, interpolation_row(x, y, mx, my, MinvT, grd))
end
