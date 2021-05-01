

function stencil1d(m::Int)

    ic = Int(ceil(m/2.0))
    M = zeros((m, m))
    for i=1:m
        for j=1:m
            M[i, j] = (i - ic)^(j-1) / factorial(j-1)
        end
    end
    return inv(M)
end


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
            M[k,:] = spline2d(i-ic, j-jc, 1.0, 1.0, mx, my)
        end
    end

    return transpose(inv(M))
end


function interpolation_row(x, y, MinvT, Nx, Ny, mx, my)

    # (i,j) indices for center point of stencil
    ics = Int(ceil(mx/2.0))
    jcs = Int(ceil(my/2.0))

    # indices for center point in grid
    icg = Int(floor(x*Nx))
    jcg = Int(floor(y*Ny))

    # values for row
    xr = x*Nx - icg
    yr = y*Ny - jcg
    spline = spline2d(xr, yr, mx, my)
    row_dat = MinvT * spline

    # collective index for first nonzero element of output
    k0 = (icg - ics + 1) + (jcg - jcs)*Nx

    # indices for nonzero elements of row
    row_js = zeros(Int, mx*my)
    for i=1:mx
        for j=1:my
            row_js[i+(j-1)*mx] = k0 + (i-1) + (j)*Nx
        end
    end

    return row_dat, row_js
end
