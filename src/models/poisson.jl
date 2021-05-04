


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


function solve_with_dirichlet(x0, f, N, a, b, dr)


    grd = Circle(N, a)
    R = outer_reflection(b, dr, grd)
    P = outer_penalization(b, dr, grd)
    L = laplacian(grd)
    fvec = f_to_grid(f, grd)
    x0vec = f_to_grid(x0, grd)

    A = L / grd.dx^2
    # A = ((I - P) * L / grd.dx^2 + P * (R + I))
    b = (I - P) * fvec

    # xa = A \ b
    # x = A \ fvec

    # xb = p_jacobi(A, x0vec, b, 1.0, 10000, 1, 2, 1e-8)
    # println(norm(xa - xb))

    x = p_jacobi(A, x0vec, fvec, 1.0, 100, 1, 2, 1e-8)

    return vec_to_mesh(x, grd)
end
