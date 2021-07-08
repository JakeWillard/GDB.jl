
function rk4_step_3d(x, y, z, fx, fy, fz, ds)

    kx1 = fx(x, y)
    ky1 = fy(x, y)
    kz1 = fz(x, y)

    kx2 = fx(x + ds * kx1 / 2, y + ds * ky1 / 2)
    ky2 = fy(x + ds * kx1 / 2, y + ds * ky1 / 2)
    kz2 = fz(x + ds * kx1 / 2, y + ds * ky1 / 2)

    kx3 = fx(x + ds * kx2 / 2, y + ds * ky2 / 2)
    ky3 = fy(x + ds * kx2 / 2, y + ds * ky2 / 2)
    kz3 = fz(x + ds * kx2 / 2, y + ds * ky2 / 2)

    kx4 = fx(x + ds * kx3, y + ds * ky3)
    ky4 = fy(x + ds * kx3, y + ds * ky3)
    kz4 = fz(x + ds * kx3, y + ds * ky3)

    x += ds * (kx1 + 2 * kx2 + 2 * kx3 + kx4) / 6
    y += ds * (ky1 + 2 * ky2 + 2 * ky3 + ky4) / 6
    z += ds * (kz1 + 2 * kz2 + 2 * kz3 + kz4) / 6

    return Float64[x, y, z]
end


function rk4_step_2d(x, y, fx, fy, ds)

    kx1 = fx(x, y)
    ky1 = fy(x, y)

    kx2 = fx(x + ds * kx1 / 2, y + ds * ky1 / 2)
    ky2 = fy(x + ds * kx1 / 2, y + ds * ky1 / 2)

    kx3 = fx(x + ds * kx2 / 2, y + ds * ky2 / 2)
    ky3 = fy(x + ds * kx2 / 2, y + ds * ky2 / 2)

    kx4 = fx(x + ds * kx3, y + ds * ky3)
    ky4 = fy(x + ds * kx3, y + ds * ky3)

    x += ds * (kx1 + 2 * kx2 + 2 * kx3 + kx4) / 6
    y += ds * (ky1 + 2 * ky2 + 2 * ky3 + ky4) / 6

    return Float64[x, y]
end


function trace_fieldline(x0, y0, bx, by, bz, ds, deltaPhi)

    x = x0
    y = y0
    z = 0
    deltaS = 0

    while z < deltaPhi
        x, y, z = rk4_step_3d(x, y, z, bx, by, bz, ds)
        deltaS += ds
    end

    return Float64[x, y, deltaS]
end



function trace_reflection(x0, y0, psi, bx, by, psi_b, ds)

    # c determines the rotation of b
    if psi(x0, y0) < psi_b
        c = 1
    else
        c = -1
    end

    fx(x, y) = c*by(x, y)
    fy(x, y) = -c*bx(x, y)

    x = x0
    y = y0
    deltaS = 0

    while c*(psi(x, y) - psi_b) <= 0

        x, y = rk4_step_2d(x, y, fx, fy, ds)
        deltaS += ds
    end

    # after reaching boundary, continue tracing distance equal to distance already traced.
    while deltaS >= 0

        x, y = rk4_step_2d(x, y, fx, fy, ds)
        deltaS -= ds
    end

    return Float64[x, y]
end
