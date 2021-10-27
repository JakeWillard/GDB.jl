


function nearest_neighbor(grd::Grid, Ns)

    _, proj_j, _ = findnz(grd.Proj)
    dK = proj_j[end] - proj_j[1]

    is = Int64[]
    js = Int64[]
    dats = Int64[]

    for j=1:grd._Ny
        for i=1:grd._Nx

            r = grd.r0 + grd.dr*[i-1, j-1]

            k_cart = i + (j-1)*grd._Nx
            if k_cart > proj_j[1]
                k_approx = Int64(ceil((k_cart - proj_j[1]) * grd.Nk / dK))
            else
                k_approx = 1
            end

            kmax = minimum([k_approx+Ns, grd.Nk])
            kmin = maximum([k_approx-Ns, 1])
            dists = [norm(r - grd.points[:,i]) for i=kmin:kmax]

            try
                mindist, _knearest = findmin(dists)
                k_nearest = [kmin:kmax...][_knearest]
                is = [is; k_cart]
                js = [js; k_nearest]
                dats = [dats; 1]
            catch
                println(k_approx)
                println(kmin)
                println(kmax)
            end
        end
    end

    return sparse(is, js, dats, grd._Nx*grd._Ny, grd.Nk)
end
