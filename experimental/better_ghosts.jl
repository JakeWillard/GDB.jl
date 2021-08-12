

struct BetterGhosts

    Proj :: SparseMatrixCSC
    Extr :: SparseMatrixCSC
    swaps :: Matrix{Float64}

end


function simple_reflection(mx, my, MinvT, barr::Barrier, grd::Grid)

    is = Int64[]
    js = Int64[]
    dat = Float64[]

    for k=1:grd.Nk

        if smoothstep(grd.points[:,k]..., 1.0, barr) < 0.5
            x, y = barr.rmap(grd.points[:,k]...)
            row_dat, row_js = interpolation_row(x, y, mx, my, MinvT, grd)
            js = [js; row_js]
            dat = [dat; row_dat]
        else
            row_dat, row_js = interpolation_row(grd.points[:,k]..., mx, my, MinvT, grd)
            js = [js; row_js]
            dat = [dat; row_dat]
        end
        is = [is; k*ones(mx*my)]
    end

    return sparse(is, js, dat, grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
end


function taxi_car_reflection(func::Function, grd::Grid)

    _, js, _ = findnz(grd.Proj)
    is, ks, _ = findnz(transpose(grd.Proj))
    jnew = Int64[]

    for k=1:grd.Nk
        x, y = grd.points[:,k]
        kc = js[k]
        i = rem(kc, grd._Nx)
        j = div(kc, grd._Ny)
        steps = 0

        while func(x,y) < 0

            fs = hcat([[func(x+grd.dx*i, y+grd.dy*j) for i=-1:1] for j=-1:1]...)
            di, dj = Tuple(findmax(fs)[2])
            di -= 2
            dj -= 2
            x += di*grd.dx
            y += dj*grd.dy
            i += di
            j += dj
            steps += 1
        end

        for _=1:steps

            fs = hcat([[func(x+grd.dx*i, y+grd.dy*j) for i=-1:1] for j=-1:1]...)
            di, dj = Tuple(findmax(fs)[2])
            di -= 2
            dj -= 2
            x += di*grd.dx
            y += dj*grd.dy
            i += di
            j += dj
        end

        knew = ks[findall(x->x==i + (j-1)*grd._Nx, is)[1]]
        jnew = [jnew; knew]
    end

    return sparse([1:grd.Nk...], jnew, ones(grd.Nk), grd.Nk, grd.Nk)
end
