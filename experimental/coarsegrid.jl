
function coarsen_grid(grd::Grid)

    _, ks, _ = findnz(grd.Proj)

    Nxc = Int64(ceil(grd._Nx/2))
    Nyc = Int64(ceil(grd._Ny/2))

    points = zeros(2, grd.Nk)
    projc_js = zeros(grd.Nk)
    Nkc = 0
    for i=1:grd.Nk
        k = ks[i]
        j = Int(ceil(k / grd._Nx))
        i = k - (j - 1) * grd._Nx
        if iseven(i) && iseven(j)
            Nkc += 1
            points[:,Nkc] = grd.points[:,k]

            ic = (i + 1)/2
            jc = (j + 1)/2
            projc_js[Nkc] = ic + (jc-1)*Nxc
        end
    end

    Projc = sparse([1:Nkc...], projc_js[1:Nkc], ones(Nkc), Nkc, Nxc*Nyc)
    return Grid(grd.r0, grd.r1, points[:,1:Nkc], Projc, grd.dr*2, Nkc, Nxc, Nyc, grd._Nbuffer, grd._nan_outside_boundaries)
end
