

function mirror_image(k, func::Function, grd::Grid)

    _, proj_js, _ = findnz(grd.Proj)

    x, y = grd.points[:,k]
    k_cart = proj_js[k]
    i = rem(k_cart, grd._Nx)
    j = div(k_cart, grd._Ny)
    steps = 0

    while func(x,y) < 0.0
        fs = hcat([[func(x+grd.dx*i, y+grd.dy*j) for i=-1:1] for j=-1:1]...)
        di, dj = Tuple(findmax(fs)[2]) .- 2
        x += di*grd.dx
        y += dj*grd.dy
        i += di
        j += dj
        steps += 1
    end

    for _=1:steps
        fs = hcat([[func(x+grd.dx*i, y+grd.dy*j) for i=-1:1] for j=-1:1]...)
        di, dj = Tuple(findmax(fs)[2]) .- 2
        x += di*grd.dx
        y += dj*grd.dy
        i += di
        j += dj
    end

    k_cart_new = i + (j-1)*grd._Nx
    knew = findall(x->x==k_cart_new, proj_js)[1]

    return knew
end


function mirror_matrix(fs::Vector{Function}, grd::Grid)

    js = DArray((grd.Nk,), workers(), length(workers())) do inds

        ks = Int64[1:grd.Nk...]
        j = Int64[]
        for i=1:length(ks)
            x, y = grd.points[:,ks[i]]
            flags = Bool[f(x,y) < 0.0 for f in fs]
            trueflag = findfirst(x->x, flags)

            f(x,y) = isnothing(trueflag) ? fs[1](x,y) : fs[trueflag](x,y)
            j = [j; mirror_image(ks[i], f, grd)]
        end

        j
    end

    return sparse([1:grd.Nk...], Array(js), ones(grd.Nk), grd.Nk, grd.Nk)
end
