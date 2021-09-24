

struct GhostData

    Proj :: SparseMatrixCSC
    R :: SparseMatrixCSC
    flip_factors :: Matrix{Int64}

end


function GhostData(M::Mirror, grd::Grid)

    chnl = RemoteChannel(()->Channel{Bool}(), 1)
    p = Progress(grd.Nk, desc="Resolving mirrors... ")
    @async while take!(chnl)
        next!(p)
    end

    map_out = pmap(1:grd.Nk) do k

        _, proj_js, _ = findnz(grd.Proj)
        k_cart = proj_js[k]
        h, sindex = distance_to_mirror(grd.points[:,k]..., M)
        flip_row = ones(Int64, size(M.verts)[2])

        if h < 0

            flip_row[sindex] = -1

            j = Int(ceil(k_cart / grd._Nx))
            i = k_cart - (j - 1) * grd._Nx
            x, y, di, dj = mirror_image(grd.points[:,k]..., grd.dx, grd.dy, M)
            i += di
            j += dj
            # @assert distance_to_mirror(x, y, M)[1] > 0
            @assert !((di == 0) && (dj == 0))
            k_cart = i + (j-1)*grd._Nx
            k = 0
        end

        # @info k
        put!(chnl, true)
        Int64[k, k_cart, flip_row...]
    end

    C = hcat(map_out...)

    proj_j = C[1,:]
    proj_j = proj_j[proj_j .!= 0]
    r_j = C[2,:]
    flip_factors = C[3:end,:]

    Mirror = sparse([1:grd.Nk...], r_j, ones(grd.Nk), grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
    Proj = sparse([1:length(proj_j)...], proj_j, ones(length(proj_j)), length(proj_j), grd.Nk)

    return GhostData(Proj, Mirror, flip_factors)
end


function flip_segments(gd::GhostData, bdry)

    R = gd.R
    fac = gd.flip_factors[bdry]
    for j=1:size(factors)[1]
        R = Diagonal(vec(fac[j,:])) * R
    end

    return GhostImage(gd.Proj, R, gd.flip_factors)
end


function extrapolate_values(x::Vector{Float64}, xb::Vector{Float64}, gd::GhostData)

    return gd.R*transpose(gc.Proj)*x + (I - gd.R)*xb
end


function require_boundary_conditions(A::SparseMatrixCSC, gd::GhostData)

    P = gd.Proj
    Pt = transpose(P)

    Anew = P*A*gd.R*Pt
    Bterm = P*A*(I - gd.R)

    return Anew, Bterm
end