

function illegal_psi_variation(sol::SOL, Nr, Nt)

    rs = LinRange(sol.fm.r0, 0.8, Nr)
    ts = LinRange(-pi, pi, Nt+1)[1:Nt]
    var = zeros(Nr, Nt)

    for i=1:Nr
        for j=1:Nt
            z0 = rs[i]*exp(im*ts[j])
            z_fwd, z_bck, _, _ = sol.fm(z0)
            p0 = sol.psi(z0)
            dp1 = sol.psi(z_fwd) - p0
            dp2 = sol.psi(z_bck) - p0
            var[j, i] = norm([dp1, dp2])
        end
    end

    return rs, ts, var
end


function trace_via_fluxmap(z0, N, fm::FluxMap)

    zs = Vector{Complex{Float64}}(undef, N)
    zs[1] = z0
    for i=2:N
        zs[i] = fm(zs[i-1])[1]
    end

    return zs
end


function bilinear_coefficients(x, y, grd::Grid)

    i = Int64(floor((x - grd.r0[1]) / grd.dr)) + 1
    j = Int64(floor((y - grd.r0[2]) / grd.dr)) + 1
    u = (x - grd.r0[1]) / grd.dr - i + 1
    v = (y - grd.r0[2]) / grd.dr - j + 1

    k0 = i + (j - 1)*grd._Nx
    ks = [k0, k0+grd._Nx, k0+1, k0+grd._Nx+1]
    C = [(1 - u)*(1 - v), (1 - u)*v, u*(1 - v), u*v]

    return ks, C
end


function parallel_derivative(fm::FluxMap, grd::Grid)

    is = Int64[]
    js = Int64[]
    dats = Float64[]

    for k=1:grd.Nk
        z0 = Complex(grd.points[:,k]...)
        z_fwd, z_bck, ds1, ds2 = fm(z0)
        ds = ds1 + ds2

        # forward part
        ks, C = bilinear_coefficients(real(z_fwd), imag(z_fwd), grd)
        is = [is; k*ones(4)]
        js = [js, ks]
        dats = [dats; C / ds]

        # backward part
        ks, C = bilinear_coefficients(real(z_bck), imag(z_bck), grd)
        is = [is; k*ones(4)]
        js = [js, ks]
        dats = [dats; C / ds]
    end

    return sparse(is, js, dats, grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
end
