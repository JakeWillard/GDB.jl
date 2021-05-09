

function Circle(N, a)

    r(x, y) = sqrt((x - 0.5)^2 + (y - 0.5)^2)
    inside(x, y) = r(x,y) <= a
    return Grid(inside, N, N, 3, 3)
end



function outer_reflection(b, dr, grd::Grid)

    pts = zeros(Float64, size(grd.points))

    for k=1:grd.Nk
        x = grd.points[1,k]
        y = grd.points[2,k]
        r = sqrt((x-0.5)^2 + (y-0.5)^2)
        theta = atan(y-0.5, x-0.5)
        if r > (b - dr)
            rp = 2*b - r
            xp = rp * cos(theta) + 0.5
            yp = rp * sin(theta) + 0.5
            pts[:,k] = [xp, yp]
        else
            pts[:,k] = [x, y]
        end
    end

    return interpolation_matrix(pts, grd)
end


function outer_penalization(b, dr, grd::Grid)

    penvec = zeros(Float64, grd.Nk)
    for k=1:grd.Nk
        x = grd.points[1,k]
        y = grd.points[2,k]
        r = sqrt((x-0.5)^2 + (y-0.5)^2)
        l = (r - b - dr)/(2*dr)
        if l < 0
            penvec[k] = 0
        elseif 0 < l < 1
            penvec[k] = 3*l^2 - 2*l^3
        else
            penvec[k] = 1
        end
    end

    return Diagonal(penvec)
end



function vorticity(w, phib, N)

    grd = Circle(N, 0.45)

    pb = f_to_grid(phib, grd)
    f = f_to_grid(w, grd)

    R = outer_reflection(0.4, 0.01, grd)
    P = outer_penalization(0.4, 0.01, grd)
    L = laplacian(grd) / grd.dx^2

    A = (I - P) * L + P * (R + I)
    b = (I - P) * f + P * (R + I) * pb

    b = transpose(A) * b
    A = transpose(A) * A

    x = A \ b
    return vec_to_mesh(x, grd)
end

function n_diff(lnn, D, N)

    grd = Circle(N, 0.45)
    f = f_to_grid(lnn, grd)

    R = outer_reflection(0.3, 0.03, grd)
    P = outer_penalization(0.3, 0.03, grd)
    L = laplacian(grd) / grd.dx^2

    M = (I - P) + P*(R - I)
    c = (I - P) * f
    c = transpose(M) * c
    M = transpose(M) * M

    f = M \ c
    return vec_to_mesh(f, grd)

    A = (I - P) * (I - D*L) + P * (R - I)
    b = (I - P) * f + 2*log(0.2)*diag(P)

    b2 = transpose(A) * b
    A2 = transpose(A) * A

    x = A2 \ b2
    return vec_to_mesh(x, grd)
end
