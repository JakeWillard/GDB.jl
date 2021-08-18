
function distance_to_segment(x, y, vert0::Vector{Float64}, vert1::Vector{Float64})

    xs, ys = Float64[x, y] - vert0
    t = begin dv = vert1 - vert0; atan(dv[2], dv[1]) end

    u = xs * cos(t) + ys * sin(t)
    v = -xs * sin(t) + ys * cos(t)

    if u < 0
        return norm([x,y] - vert0) * sign(v)
    elseif 0 <= u < norm(vert1 - vert0)
        return v
    else
        return norm([x,y] - vert1) * sign(v)
    end
end


function mirror_image(x, y, dx, dy, vert0::Vector{Float64}, vert1::Vector{Float64})

    i = 0
    j = 0
    dist = 0
    steps = 0

    while distance_to_segment(x, y, vert0, vert1) < 0
        hs = hcat([[distance_to_segment(x+ii*dx, y+jj*dy, vert0, vert1) for ii=-1:1] for jj=-1:1]...)
        di, dj = Tuple(findmax(hs)[2]) .- 2
        x += di*grd.dx
        y += dj*grd.dy
        i += di
        j += dj
        dist += sqrt(di^2 + dj^2)
        steps += di^2 + dj^2
    end

    while (dist > 0) && (steps > 0)
        hs = hcat([[distance_to_segment(x+ii*dx, y+jj*dy, vert0, vert1) for ii=-1:1] for jj=-1:1]...)
        di, dj = Tuple(findmax(hs)[2]) .- 2
        x += di*grd.dx
        y += dj*grd.dy
        i += di
        j += dj
        dist -= sqrt(di^2 + dj^2)
        steps -= 1
    end

    return x, y, i, j
end


struct MirrorSegments
    verts :: Array{Float64, 3}
    orientations :: Vector{Int64}
end


@recipe function f(ms::MirrorSegments)

    Ns = size(ms.verts)[3]
    x_vectors = Vector{Float64}[]
    y_vectors = Vector{Float64}[]
    xvec = Float64[]
    yvec = Float64[]

    for i=1:Ns-1

        r0 = ms.verts[:,1,i]
        r1 = ms.verts[:,2,i]
        r2 = ms.verts[:,1,i+1]
        xvec = [xvec; r0[1]]
        yvec = [yvec; r0[2]]

        if all(r1 .!= r2)
            xvec = [xvec; r1[1]]
            yvec = [yvec; r1[2]]
            append!(x_vectors, [xvec])
            append!(y_vectors, [yvec])
            xvec = Float64[]
            yvec = Float64[]
        end
    end

    # add the last segment
    xvec = [xvec; ms.verts[1,:,end]]
    yvec = [yvec; ms.verts[2,:,end]]
    append!(x_vectors, [xvec])
    append!(y_vectors, [yvec])

    Nc = length(x_vectors)
    sizes = [length(v) for v in x_vectors]

    x = fill(NaN, (Ns, Nc))
    y = fill(NaN, (Ns, Nc))

    for i=1:Nc
        x[1:sizes[i],i] = x_vectors[i][:]
        y[1:sizes[i],i] = y_vectors[i][:]
    end

    return x, y
end


function distance_to_mirror(x, y, ms::MirrorSegments)

    Ns = size(ms.verts)[3]
    hs = Float64[ms.orientations[k] * distance_to_segment(x, y, ms.verts[:,1,k], ms.verts[:,2,k]) for k=1:Ns]
    _, closest_segment = findmin(hs .^ 2)

    return hs[closest_segment], closest_segment
end


function mirror_image(x, y, dx, dy, ms::MirrorSegments)

    # find closest segment
    _, k = distance_to_mirror(x, y, ms)
    vert0 = ms.verts[:,1,k]
    vert1 = ms.verts[:,2,k]

    return mirror_image(x, y, dx, dy, vert0, vert1)
end


struct GhostImage

    Proj :: SparseMatrixCSC
    Img :: SparseMatrixCSC
    flip_factors :: Matrix{Int64}

end


function GhostImage(ms::MirrorSegments, grd::Grid)

    chnl = RemoteChannel(()->Channel{Bool}(), 1)
    p = Progress(grd.Nk, desc="Resolving mirrors... ")
    @async while take!(chnl)
        next!(p)
    end

    map_out = pmap(1:grd.Nk) do k

        _, proj_js, _ = findnz(grd.Proj)
        k_cart = proj_js[k]
        h, sindex = distance_to_mirror(grd.points[:,k]..., ms)
        flip_row = ones(Int64, size(ms.verts)[3])

        if h < 0

            flip_row[sindex] = -1

            j = Int(ceil(k_cart / grd._Nx))
            i = k_cart - (j - 1) * grd._Nx
            x, y, di, dj = mirror_image(grd.points[:,k]..., grd.dx, grd.dy, ms)
            i += di
            j += dj
            k_cart = i + (j-1)*grd._Nx
            k = 0
        end

        put!(chnl, true)
        Int64[k, k_cart, flip_row...]
    end

    M = hcat(map_out...)

    proj_j = M[1,:]
    proj_j = proj_j[proj_j .!= 0]
    img_j = M[2,:]
    flip_factors = M[3:end,:]

    Img = sparse([1:grd.Nk...], img_j, ones(grd.Nk), grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
    Proj = sparse([1:length(proj_j)...], proj_j, ones(length(proj_j)), length(proj_j), grd.Nk)

    return GhostImage(Proj, Img, flip_factors)
end


function flip_segments(gi::GhostImage, ranges...)

    Img = gi.Img

    for rng in ranges
        factors = gi.flip_factors[rng,:]
        for j=1:size(factors)[1]
            Img = Diagonal(vec(factors[j,:])) * Img
        end
    end

    return GhostImage(gi.Proj, Img, gi.flip_factors)
end
