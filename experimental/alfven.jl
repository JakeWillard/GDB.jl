



function bilinear_coefficients(x, y, grd::Grid)

    i = Int64(floor((x - grd.r0[1]) / grd.dr)) + 1 + grd._Nbuffer
    j = Int64(floor((y - grd.r0[2]) / grd.dr)) + 1 + grd._Nbuffer
    u = rem(x - grd.r0[1], grd.dr) / grd.dr
    v = rem(y - grd.r0[2], grd.dr) / grd.dr

    k0 = i + (j - 1)*grd._Nx
    ks = Int64[k0, k0+grd._Nx, k0+1, k0+grd._Nx+1]
    C = [(1 - u)*(1 - v), (1 - u)*v, u*(1 - v), u*v]

    return ks, C
end


function parallel_laplacian(ds, Nz::Int64, sol::SOL, dm::DMap, grd::Grid)

    is = Int64[]
    js = Int64[]
    dats = Float64[]
    dPhi = 2*pi / Nz

    p = Progress(grd.Nk)

    for k=1:grd.Nk
        z0 = dm(grd.points[:,k]...)
        z_fwd, ds_fwd = fluxmap(z0, dPhi, ds, sol)
        z_bck, ds_bck = fluxmap(z0, dPhi, -ds, sol)
        ds = (ds_fwd + ds_bck)/2

        ks, C = bilinear_coefficients(dm(z0)..., grd)
        is = [is; k*ones(Int64, 4)]
        js = [js; ks]
        dats = [-2*dats; C / ds^2]

        ks, C = bilinear_coefficients(dm(z_fwd)..., grd)
        is = [is; k*ones(Int64, 4)]
        js = [js; ks]
        dats = [dats; C / ds^2]

        ks, C = bilinear_coefficients(dm(z_bck)..., grd)
        is = [is; k*ones(Int64, 4)]
        js = [js; ks]
        dats = [dats; C / ds^2]

        next!(p)
    end
    # eturn is, js, dats
    return is, js, dats


    return sparse(is, js, dats, grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
end
