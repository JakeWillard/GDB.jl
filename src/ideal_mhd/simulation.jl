


"""

INITIAL CONDITIONS

LogDensity: 0.15 everywhere, but 0.3 on x_flux surface and 0 on y_flux surface
Vorticity: 0 everywhere
Flux: 0 everywhere
Stream: kx*ky*x*y (actual initial value will always be determined by vorticity equation, but we
                   initialize this way since that also sets the boundary value.)
Current: 0 everywhere (consistent with flux)

"""


function LogDensity(stp::Setup)

    val = 0.15*ones(Float64, stp.grd.Nk)
    val = stp.P1*val[:] + 0.3*diag(I - stp.P1)
    val = stp.P2*val[:]
    return LogDensity(val)
end


Vorticity(stp::Setup) = Vorticity(zeros(Float64, stp.grd.Nk))
Flux(stp::Setup) = Flux(zeros(Float64, stp.grd.Nk))
Current(stp::Setup) = Current(zeros(Float64, stp.grd.Nk))


function Stream(kx, ky, stp)

    val = f_to_grid((x,y) -> kx*ky*x*y, stp.grd)
    return Stream(val)
end


function initialize_output_file(path, N, kx, ky, stp)

    lnn = LogDensity(stp)
    w = Vorticity(stp)
    psi = Flux(stp)
    phi = Stream(kx, ky, stp)
    j = Current(stp)

    fid = h5open(path, "w")

    fid["LogDensity"] = zeros(Float64, (stp.grd.Nk, N))
    fid["Flux"] = zeros(Float64, (stp.grd.Nk, N))
    fid["Stream"] = zeros(Float64, (stp.grd.Nk, N))
    fid["Vorticity"] = zeros(Float64, (stp.grd.Nk, N))
    fid["Current"] = zeros(Float64, (stp.grd.Nk, N))

    fid["LogDensity"][:,1] = lnn[:,2]
    fid["Fux"][:,t] = psi[:,1]
    fid["Stream"][:,t] = phi[:,1]
    fid["Vorticity"][:,t] = w[:,1]
    fid["Current"][:,t] = j[:,1]
    fid["t"] = 2

    close(fid)
end


function simulate(init_path, stp_path, dt, eta, Cdiff, n)

    fid_init = h5open(init_path, "r")
    lnn = LogDensity(fid["LogDensity"][:,1])
    psi = Flux(fid["Flux"][:,1])
    phi = Stream(fid["Stream"][:,1])
    w = Vorticity(fid["Vorticity"][:,1])
    j = Current(fid["Current"][:,1])
    N = size(fid["LogDensity"])[2]
    close(fid_init)

    fid_stp = h5open(stp_path, "r")
    stp = load_setup(fid_stp)
    close(fid_stp)

    for dummy=1:N
        integrate_then_save!(init_path, n, lnn, w, psi, phi, j, phi, dt, eta, Cdiff, stp)
    end
end
