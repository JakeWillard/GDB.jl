

abstract type LogDensity <: Physical end
abstract type Vorticity <: Physical end
abstract type Potential <: Physical end


function Annalus(N, a, b)

    r(x, y) = sqrt((x - 0.5)^2 + (y - 0.5)^2)
    inside(x, y) = a <= r(x,y) <= b
    return Grid(inside, N, N, 3, 3)
end


function inner_reflection(a, dr, grd::Grid)

    pts = zeros(Float64, size(grd.points))

    for k=1:grd.Nk
        x = grd.points[1,k]
        y = grd.points[2,k]
        r = sqrt((x-0.5)^2 + (y-0.5)^2)
        theta = atan(y-0.5, x-0.5)
        if r < (a + dr)
            rp = 2*a - r
            xp = rp * cos(theta) + 0.5
            yp = rp * sin(theta) + 0.5
            pts[:,k] = [xp, yp]
        else
            pts[:,k] = [x, y]
        end
    end

    return interpolation_matrix(pts, grd)
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


# NOTE: uses "smoothstep" function found on wikipedia
function inner_penalization(a, dr, grd::Grid)

    penvec = zeros(Float64, grd.Nk)
    for k=1:grd.Nk
        x = grd.points[1,k]
        y = grd.points[2,k]
        r = sqrt((x-0.5)^2 + (y-0.5)^2)
        l = (a + dr - r)/(2*dr)
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


function partial_t!(w::Variable{Vorticity}, phi::Variable{Potential}, R1, R2, P1, P2)

    w.t = phi.x[:,:] .* w.y[:,:] - phi.y[:,:] .* w.x[:,:]
    # w.t = (I - P1) * w.t + P1 * w.cval
    # w.t = (I - P2) * w.t + P2 * w.cval
    w.t = (I - P1) * w.t - P1 * (I + R1) * (10*w.cval)
    w.t = (I - P2) * w.t - P2 * (I + R2) * (10*w.cval)
end


function partial_t!(lnn::Variable{LogDensity}, phi::Variable{Potential}, R1, R2, P1, P2)

    lnn.t = phi.x[:,:] .* lnn.y[:,:] - phi.y[:,:] .* lnn.x[:,:]
    # lnn.t = (I - P1) * lnn.t[:,:] - P1 * (I + R1) * (lnn.cval[:,:] .- log(0.2))
    # lnn.t = (I - P2) * lnn.t[:,:] - P2 * (I - R2) * lnn.cval[:,:]
end


function apply_diffusion!(w::Variable{Vorticity}, D, L, R1, R2, P1, P2)

    A = I - D*L
    A = (I - P1)*A + P1*(R1 + I)
    A = (I - P2)*A + P2*(R2 + I)
    # A += P1*(R1 + I) - P1*A
    # A += P2*(R2 + I) - P2*A
    x0 = w.fval[:,1]
    b = (I - P1) * x0
    b -= P2 * b[:]

    w.fval[:,1] = p_jacobi(A, x0, b, 0.2, 100, 1, 2, 1e-8)
end


function apply_diffusion!(lnn::Variable{LogDensity}, D, L, R1, R2, P1, P2)

    A = I - D*L
    A += 0.5*P1*(R1 + I) - P1*A
    A += P2*(R2 - I) - P2*A
    x0 = exp.(lnn.fval[:,1])
    b = (I - P1) * x0[:] + P1 * (log(0.2)*ones(Float64, size(x0)))
    b = (I - P2) * b[:]

    lnn.fval[:,1] = p_jacobi(A, x0, b, 0.2, 100, 1, 2, 1e-8)
end


function solve_vorticity_eqn!(phi::Variable{Potential}, lnn::Variable{LogDensity}, w::Variable{Vorticity}, phi0::Matrix{Float64}, L, P1, P2)

    A = (I - P1) * L + P1
    A += P2 - P1*A
    b = w.fval[:,1] # .* exp.(-lnn.fval[:,1])
    b[:] = (I - P1) * b[:] + P1 * phi0[:,1]
    b[:] = (I - P2) * b[:] + P2 * phi0[:,2]

    phi.cval[:,1] = p_jacobi(A, phi.cval[:,1], b, 0.2, 100, 1, 2, 1e-8)
end



function write_data(path, lnn, w, phi, t)

    fid = h5open(path, "r+")
    fid["LogDensity"][:,t] = lnn.cval[:,1]
    fid["Vorticity"][:,t] = w.cval[:,1]
    fid["Potential"][:,t] = phi.cval[:,1]
    close(fid)
end



function timestep!(lnn, w, phi, grd, phi0, dt, D, Dx, Dy, L, R1, R2, P1, P2)

    phi.x = Dx * phi.cval / grd.dx
    phi.y = Dy * phi.cval / grd.dy
    lnn.x = Dx * lnn.cval / grd.dx
    lnn.y = Dy * lnn.cval / grd.dy
    w.x = Dx * w.cval / grd.dx
    w.y = Dy * w.cval / grd.dy
    partial_t!(lnn, phi, R1, R2, P1, P2)
    partial_t!(w, phi, R1, R2, P1, P2)

    lnn.fval = 0.5*(lnn.pval[:,:] + lnn.cval[:,:]) + dt*lnn.t[:,:]
    w.fval = 0.5*(w.pval[:,:] + w.cval[:,:]) + dt*w.t[:,:]

    # apply_diffusion!(lnn, D, L, R1, R2, P1, P2)
    apply_diffusion!(w, D, L, R1, R2, P1, P2)
    # solve_vorticity_eqn!(phi, lnn, w, phi0, L, P1, P2)
    lnn.pval = lnn.cval[:,:]
    lnn.cval = lnn.fval[:,:]
    w.pval = w.cval[:,:]
    w.cval = w.fval[:,:]

    phi.x = Dx * phi.cval[:,:] / grd.dx
    phi.y = Dy * phi.cval[:,:] / grd.dy
    lnn.x = Dx * lnn.cval[:,:] / grd.dx
    lnn.y = Dy * lnn.cval[:,:] / grd.dy
    w.x = Dx * w.cval[:,:] / grd.dx
    w.y = Dy * w.cval[:,:] / grd.dy
    partial_t!(lnn, phi, R1, R2, P1, P2)
    partial_t!(w, phi, R1, R2, P1, P2)

    lnn.fval = lnn.pval[:,:] + dt*lnn.t[:,:]
    w.fval = w.pval[:,:] + dt*w.t[:,:]

    # apply_diffusion!(lnn, D, L, R1, R2, P1, P2)
    apply_diffusion!(w, D, L, R1, R2, P1, P2)
    # solve_vorticity_eqn!(phi, lnn, w, phi0, L, P1, P2)
    lnn.pval = lnn.cval[:,:]
    lnn.cval = lnn.fval[:,:]
    w.pval = w.cval[:,:]
    w.cval = w.fval[:,:]
end


function batch_integrate(t, path, lnn, w, phi, args...; nt=10)

    for i=1:nt
        timestep!(lnn, w, phi, args...)
    end

    write_data(path, lnn, w, phi, t)
end


function simulation(path, N, Nt)

    dt = 0.001
    D = 0.001

    grd = Annalus(N, 0.1, 0.45)

    Dx = x_derivative(1, grd)
    Dy = y_derivative(1, grd)
    L = laplacian(grd)

    R1 = inner_reflection(0.12, 0.01, grd)
    R2 = outer_reflection(0.42, 0.01, grd)
    P1 = inner_penalization(0.12, 0.01, grd)
    P2 = outer_penalization(0.42, 0.01, grd)

    phi0 = zeros(Float64, (grd.Nk, 2))
    phi0[:,1] = f_to_grid((x,y) -> y^2, grd)

    lnn = Variable{LogDensity}((x,y) -> log(0.2), 1, grd)
    w = Variable{Vorticity}((x,y) -> sin(3*pi*x)*sin(3*pi*y), 1, grd)
    phi = Variable{Potential}((x,y) -> 3*sqrt((x-0.5)^2 + (y-0.5)^2), 1, grd)

    fid = h5open(path, "w")

    save_grid(fid, grd)

    fid["LogDensity"] = zeros(Float64, (grd.Nk, Nt))
    fid["Vorticity"] = zeros(Float64, (grd.Nk, Nt))
    fid["Potential"] = zeros(Float64, (grd.Nk, Nt))

    fid["LogDensity"][:,1] = lnn.cval[:,1]
    fid["Vorticity"][:,1] = w.cval[:,1]
    fid["Potential"][:,1] = phi.cval[:,1]

    close(fid)

    for t=2:Nt
        batch_integrate(t, path, lnn, w, phi, grd, phi0, dt, D, Dx, Dy, L, R1, R2, P1, P2, nt=2)
        println(t)
    end
end




# function setup_simulation(path, lnn0, w0, phi0, N, a, b, dr, dt, D)
#
#     grd = Annalus(N, a, b)
#
#     Dx = x_derivative(1, grd)
#     Dy = y_derivative(1, grd)
#     L = laplacian(grd)
#
#     R1 = inner_reflection(a, dr, grd)
#     R2 = outer_reflection(b, dr, grd)
#     P1 = inner_penalization(a, dr, grd)
#     P2 = outer_penalization(b, dr, grd)
#
#     lnn = Variable{LogDensity}(lnn0, 1, grd)
#     w = Variable{Vorticity}(w0, 1, grd)
#     phi = Variable{Potential}(phi0, 1, grd)
#
#     fid = h5open(path, "w")
#
#     gDx = create_group(fid, "Dx")
#     gDy = create_group(fid, "Dy")
#     gL = create_group(fid, "L")
#     gPinv = create_group(fid, "Pinv")
#     gP = create_group(fid, "P")
#     gR1 = create_group(fid, "R1")
#     gR2 = create_group(fid, "R2")
#     gP1 = create_group(fid, "P1")
#     gP2 = create_group(fid, "P2")
#
#     out = create_group(fid, "output")
#     misc = create_group(fid, "misc")
#
#     save_matrix!(gDx, Dx)
#     save_matrix!(gDy, Dy)
#     save_matrix!(gL, L)
#     save_matrix!(gR1, R1)
#     save_matrix!(gR2, R2)
#     save_matrix!(gP1, P1)
#     save_matrix!(gP2, P2)
#     save_matrix!(gPinv, grd.inverse_projection)
#     save_matrix!(gP, grd.projection)
#
#     misc["points"] = grd.points[:,:]
#     misc["params"] = Float64[grd.dx, grd.dy, dt, a, b, D]
#
#     out["LogDensity"] = zeros(Float64, (grd.Nk, 1, 2))
#     out["Potential"] = zeros(Float64, (grd.Nk, 1, 2))
#     out["Vorticity"] = zeros(Float64, (grd.Nk, 1, 2))
#
#     out["LogDensity"][:,:,1] = lnn.pval[:,:]
#     out["LogDensity"][:,:,2] = lnn.cval[:,:]
#     out["Vorticity"][:,:,1] = w.pval[:,:]
#     out["Vorticity"][:,:,2] = w.cval[:,:]
#     out["Potential"][:,:,1] = phi.pval[:,:]
#     out["Potential"][:,:,2] = phi.cval[:,:]
#
#     close(fid)
# end
#
#
# function continue_simulation(path, Nt)
#
#     fid = h5open(path, "r+")
#
#     Dx = read_matrix(fid["Dx"])
#     Dy = read_matrix(fid["Dy"])
#     L = read_matrix(fid["L"])
#     R1 = read_matrix(fid["R1"])
#     R2 = read_matrix(fid["R2"])
#     P1 = read_matrix(fid["P1"])
#     P2 = read_matrix(fid["P2"])
#
#     Nk, Nz, Nt0 = size(fid["output/LogDensity"][:,:,:])
#     lnn_dset = zeros(Float64, (Nk, Nz, Nt0 + Nt))
#     w_dset = zeros(Float64, (Nk, Nz, Nt0 + Nt))
#     phi_dset = zeros(Float64, (Nk, Nz, Nt0 + Nt))
#
#     lnn_dset[:,:,1:Nt0] = fid["output/LogDensity"][:,:,:]
#     w_dset[:,:,1:Nt0] = fid["output/Vorticity"][:,:,:]
#     phi_dset[:,:,1:Nt0] = fid["output/Potential"][:,:,:]
#
#     fid["output/LogDensity"] = lnn_dset[:,:,:]
#     fid["output/Vorticity"] = w_dset[:,:,:]
#     fid["output/Potential"] = phi_dset[:,:,:]
#
#     lnn = Variable{LogDensity}(lnn_dset[:,:,Nt0-1], lnn_dset[:,:,Nt0])
#
#
#
#
#
# end
