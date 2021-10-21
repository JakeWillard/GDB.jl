

struct Grid

    r0 :: Vector{Float64}
    r1 :: Vector{Float64}
    points :: Matrix{Float64}
    Proj :: SparseMatrixCSC
    dr :: Float64
    Nk :: Int64
    _Nx :: Int64
    _Ny :: Int64
    _Nbuffer :: Int64
    _nan_outside_boundaries :: Matrix{Float64}

end

function Grid(M::Mirror, h::Float64, r0::Vector{Float64}, r1::Vector{Float64}, p::Int64; Nbuffer=100)

    deltaX = r1[1] - r0[1]
    deltaY = r1[2] - r0[2]

    if deltaX < deltaY
        Nx = 2^p
        dr = deltaX / Nx
        Ny = Int64(div(deltaY, dr))
    else
        Ny = 2^p
        dr = deltaY / Ny
        Nx = Int64(div(deltaX, dr))
    end

    # recalculate r1
    r1 = r0 + dr*[Nx, Ny]

    _Nx = Nx + 2*Nbuffer
    _Ny = Ny + 2*Nbuffer

    points = zeros(Float64, (2, Nx*Ny))
    proj_rows = Int32[i for i=1:Nx*Ny]
    proj_cols = zeros(Int32, Nx*Ny)
    proj_vals = ones(Int32, Nx*Ny)
    _nan_outside_boundaries = fill(NaN, (Nx, Ny))
    k = 0

    for j=1:Ny
        for i=1:Nx

            x = r0[1] + (i-1)*dr
            y = r0[2] + (j-1)*dr
            if distance_to_mirror(x, y, M)[1] > -h
                k += 1
                points[:,k] = [x, y]
                proj_cols[k] = i + _Nx*(j-1+Nbuffer) + Nbuffer
                _nan_outside_boundaries[i,j] = 1.0
            end

        end
    end

    Proj = sparse(proj_rows[1:k], proj_cols[1:k], proj_vals[1:k], k, _Nx*_Ny)
    return Grid(r0, r1, points[:,1:k], Proj, dr, k, _Nx, _Ny, Nbuffer, _nan_outside_boundaries)
end


function coarse_grid(grd::Grid)

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


@recipe function f(v::Vector{Float64}, grd::Grid)

    vp_cart = transpose(grd.Proj) * v
    V = reshape(vp_cart, grd._Nx, grd._Ny)
    Vplot = transpose(V[grd._Nbuffer+1:grd._Nx-grd._Nbuffer, grd._Nbuffer+1:grd._Ny-grd._Nbuffer] .* grd._nan_outside_boundaries)

    Nx = grd._Nx - 2*grd._Nbuffer
    Ny = grd._Ny - 2*grd._Nbuffer
    x = LinRange(grd.r0[1], grd.r1[1], Nx)
    y = LinRange(grd.r0[2], grd.r1[2], Ny)

    xlims --> (x[1], x[end])
    ylims --> (y[1], y[end])
    size --> (1000, 1000)
    aspect_ratio --> :equal

    x, y, Vplot
end


@recipe function f(v::Vector{Float64}, grd::Grid, z::Int64)

    vp = reshape(v, grd.Nk, :)[:,z]
    vp_cart = transpose(grd.Proj) * vp
    V = reshape(vp_cart, grd._Nx, grd._Ny)
    Vplot = transpose(V[grd._Nbuffer+1:grd._Nx-grd._Nbuffer, grd._Nbuffer+1:grd._Ny-grd._Nbuffer] .* grd._nan_outside_boundaries)

    Nx = grd._Nx - 2*grd._Nbuffer
    Ny = grd._Ny - 2*grd._Nbuffer
    x = LinRange(grd.r0[1], grd.r1[1], Nx)
    y = LinRange(grd.r0[2], grd.r1[2], Ny)

    xlims --> (x[1], x[end])
    ylims --> (y[1], y[end])
    size --> (1000, 1000)
    aspect_ratio --> :equal

    x, y, Vplot
end


function function_to_grid(f::Function, grd::Grid)

    return Float64[f(grd.points[:,k]...) for k=1:grd.Nk]
end


function operator_to_grid(O::Function, grd::Grid)

    return grd.Proj * O() * transpose(grd.Proj)
end
