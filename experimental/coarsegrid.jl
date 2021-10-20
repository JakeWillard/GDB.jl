
function coarsen_grid(grd::Grid)

    _, ks, _ = findnz(grd.Proj)

    _Nxc = Int64(ceil(grd._Nx/2))
    _Nyc = Int64(ceil(grd._Ny/2))
    Nbuffer = Int64(ceil(grd._Nbuffer/2))

    dr = 2*grd.dr
    if iseven(grd._Nbuffer)
        r0 = grd.r0
        nan_outside_boundaries = grd._nan_outside_boundaries[1:2:end,1:2:end]
    else
        r0 = grd.r0 + [1,1]*grd.dr
        nan_outside_boundaries = grd._nan_outside_boundaries[2:2:end,2:2:end]
    end
    Nxc, Nyc = size(nan_outside_boundaries)
    r1 = r0 + dr*[Nxc, Nyc]

    points = zeros(2, grd.Nk)
    projc_js = zeros(grd.Nk)
    Nkc = 0
    for k=1:grd.Nk
        kc = ks[k]
        j = Int(ceil(kc / grd._Nx))
        i = kc - (j - 1) * grd._Nx
        if isodd(i) && isodd(j)
            Nkc += 1
            points[:,Nkc] = grd.points[:,k]

            ic = (i + 1)/2
            jc = (j + 1)/2
            projc_js[Nkc] = ic + (jc-1)*_Nxc
        end
    end

    Projc = sparse([1:Nkc...], projc_js[1:Nkc], ones(Nkc), Nkc, _Nxc*_Nyc)
    return Grid(r0, r1, points[:,1:Nkc], Projc, dr, Nkc, _Nxc, _Nyc, Nbuffer, nan_outside_boundaries)
end
