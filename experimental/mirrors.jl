

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


struct BetterGhostConditions

    Proj :: SparseMatrixCSC
    Mirror :: SparseMatrixCSC
    swaps :: Matrix{Int64}

end


function BetterGhostConditions(bars::Vector{Barrier}, grd::Grid)

    M = DArray((grd.Nk,), workers(), (length(workers()), length(bars)+2)) do inds

        ks = Int64[1:grd.Nk...]
        swaps = ones(Int64, (grd.Nk, length(bars)))
        j_mir = Int64[]
        j_proj = Int64[]

        for i=1:length(ks)
            x, y = grd.points[:,ks[i]]

            flags = Bool[f(x,y) < 0.0 for f in fs]
            trueflag = findfirst(x->x, flags)
            if isnothing(trueflag)
                j_mir = [j_mir; ks[i]]
                j_proj = [j_proj; ks[i]]
            else
                swaps[i,trueflag] = -1
                j_mir = [j_mir; mirror_image(ks[i], (x,y) -> bars[trueflag].func(x,y) - bars[trueflag].val, grd)]
            end

        end

        hcat(j_mir, j_proj, swaps)
    end

    j_mir = Array(M[:,1])
    j_proj = Array(M[:,2])
    swaps = Array(M[:,3:end])

    Mirror = sparse([1:grd.Nk...], j_mir, ones(grd.Nk), grd.Nk, grd.Nk)
    Proj = sparse([1:length(j_proj)...], j_proj, ones(length(j_proj)), length(j_proj), grd.Nk)

    return BetterGhostConditions(Proj, Mirror, swaps)
end


function swap_sign(gc::BetterGhostConditions, inds...)

    Mirror = gc.Mirror

    for i in inds
        Mirror = Diagonal(gc.swaps[:,i]) * Mirror
    end

    return BetterGhostConditions(gc.Proj, Mirror, gc.swaps)
end


function extrapolate_ghosts(x::Vector{Float64}, xb::Vector{Float64}, gc::BetterGhostConditions)

    return gc.Mirror*transpose(gc.Proj)*x + (I - gc.Mirror)*xb
end


function require_boundary_conditions(A::SparseMatrixCSC, gc::BetterGhostConditions)

    P = gc.Proj
    Pt = transpose(P)

    Anew = P*A*gc.Mirror*Pt
    Bterm = P*A*(I - gc.Mirror)

    return Anew, Bterm
end
