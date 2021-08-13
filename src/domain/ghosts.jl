

function mirror_image(k, b::Barrier, grd::Grid)

    _, proj_js, _ = findnz(grd.Proj)

    x, y = grd.points[:,k]
    k_cart = proj_js[k]
    j = Int(ceil(k_cart / grd._Nx))
    i = k_cart - (j - 1) * grd._Nx
    dist = 0
    maxsteps = 0

    while b.orientation*(b.func(x,y) - b.val) < 0
        fs = hcat([[b.orientation*(b.func(x+grd.dx*ii, y+grd.dy*jj) - b.val) for ii=-1:1] for jj=-1:1]...)
        di, dj = Tuple(findmax(fs)[2]) .- 2
        x += di*grd.dx
        y += dj*grd.dy
        i += di
        j += dj
        dist += sqrt(di^2 + dj^2)
        maxsteps += di^2 + dj^2
    end

    steps = 0
    while (dist > 0) && (steps < maxsteps)
        fs = hcat([[b.orientation*(b.func(x+grd.dx*ii, y+grd.dy*jj) - b.val) for ii=-1:1] for jj=-1:1]...)
        di, dj = Tuple(findmax(fs)[2]) .- 2
        x += di*grd.dx
        y += dj*grd.dy
        i += di
        j += dj
        dist -= sqrt(di^2 + dj^2)
        steps += 1
    end

    k_cart_new = i + (j-1)*grd._Nx
    # findk = findall(x->x==k_cart_new, proj_js)
    # knew = isempty(findk) ? 0 : findk[1]
    # knew = findk[1]

    return k_cart_new
end


struct GhostConditions

    Proj :: SparseMatrixCSC
    Mirror :: SparseMatrixCSC
    swaps :: Matrix{Int64}

end


function GhostConditions(bars::Vector{Barrier}, grd::Grid)

    M = DArray((grd.Nk,length(bars)+2), workers(), (length(workers()), 1)) do inds

        ks = Int64[inds[1]...]
        swaps = ones(Int64, (length(inds[1]), length(bars)))
        j_mir = Int64[]
        j_proj = Int64[]
        _, proj_js, _ = findnz(grd.Proj)

        @showprogress "Computing Mirrors..." for i=1:length(ks)
            x, y = grd.points[:,ks[i]]

            flags = Bool[b.orientation*(b.func(x,y)-b.val) < 0 for b in bars]
            trueflag = findfirst(x->x, flags)
            if isnothing(trueflag)
                j_mir = [j_mir; proj_js[ks[i]]]
                j_proj = [j_proj; ks[i]]
            else
                swaps[i,trueflag] = -1
                j_mir = [j_mir; mirror_image(ks[i], bars[trueflag], grd)]
                j_proj = [j_proj; 0]
            end

            # @info "$i/$(length(ks))"

        end

        hcat(j_mir, j_proj, swaps)
    end

    j_mir = Array(M[:,1])
    j_proj = Array(M[:,2])
    j_proj = j_proj[j_proj .!= 0]
    swaps = Array(M[:,3:end])

    Mirror = sparse([1:grd.Nk...], j_mir, ones(grd.Nk), grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
    Proj = sparse([1:length(j_proj)...], j_proj, ones(length(j_proj)), length(j_proj), grd.Nk)

    return GhostConditions(Proj, Mirror, swaps)
end


function swap_sign(gc::GhostConditions, inds...)

    Mirror = gc.Mirror

    for i in inds
        Mirror = Diagonal(gc.swaps[:,i]) * Mirror
    end

    return GhostConditions(gc.Proj, Mirror, gc.swaps)
end


function extrapolate_ghosts(x::Vector{Float64}, xb::Vector{Float64}, gc::GhostConditions)

    return gc.Mirror*transpose(gc.Proj)*x + (I - gc.Mirror)*xb
end


function require_boundary_conditions(A::SparseMatrixCSC, gc::GhostConditions)

    P = gc.Proj
    Pt = transpose(P)

    Anew = P*A*gc.Mirror*Pt
    Bterm = P*A*(I - gc.Mirror)

    return Anew, Bterm
end
