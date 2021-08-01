
struct Grid

    r0 :: Vector{Float64}
    points :: Matrix{Float64}
    Proj :: SparseMatrixCSC
    dx :: Float64
    dy :: Float64
    Nk :: Int64
    Nz :: Int64
    _Nx :: Int64
    _Ny :: Int64
    _nan_outside_boundaries :: Matrix{Float64}

end


# TODO: parallelize this calculation
function Grid(is_inside::Function, r0::Vector{Float64}, r1::Vector{Float64}, Nx, Ny)

    deltaX = r1[1] - r0[1]
    deltaY = r1[2] - r0[2]
    dx = deltaX / Nx
    dy = deltaY / Ny

    points = zeros(Float64, (2, Nx*Ny))
    proj_rows = Int32[i for i=1:Nx*Ny]
    proj_cols = zeros(Int32, Nx*Ny)
    proj_vals = ones(Int32, Nx*Ny)
    _nan_outside_boundaries = fill(NaN, (Nx, Ny))
    k = 0

    for j=1:Ny
        for i=1:Nx

            x = r0[1] + (i-1)*dx
            y = r0[2] + (j-1)*dx
            if is_inside(x, y)
                k += 1
                points[:,k] = [x, y]
                proj_cols[k] = i + Nx*(j-1)
                _nan_outside_boundaries[i,j] = 1.0
            end

        end
    end

    Proj = sparse(proj_rows[1:k], proj_cols[1:k], proj_vals[1:k], k, Nx*Ny)
    return Grid(r0, points[:,1:k], Proj, dx, dy, Nx, Ny, _nan_outside_boundaries)
end



function vec_to_mesh(vec, grd::Grid)


    vals = transpose(grd.Proj) * vec

    Nx = grd._Nx
    Ny = grd._Ny
    Nz = grd.Nz
    out = zeros(Float64, (Nx, Ny, Nz))
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                l = (i + (j-1)*Nx) + (k-1)*Nx*Ny
                out[i,j,k] = vals[l] * grd._nan_outside_boundaries[i,j]
            end
        end
    end

    return out
end


function f_to_grid(f::Function, grd::Grid)

    vec = Float64[f(grd.points[:,k]...) for k=1:grd.Nk]

    return vcat([vec for _=1:grd.Nz]...)
end


function grid_map(f::Function, outsize::Tuple{Int64, Int64}, grd::Grid)

    n, m = outsize
    output = SharedArray{Float64}((n*grd.Nk,m))

    # @distributed unfortunately only works if nprocs > 1
    if nprocs() > 1
        @distributed for i=1:grd.Nk
            k = 1 + n*(i-1)
            output[k:k+n-1,:] = f(grd.points[:,i]..., grd)
        end
    else
        for i=1:grd.Nk
            k = 1 + n*(i-1)
            output[k:k+n-1,:] = f(grd.points[:,i]..., grd)
        end
    end

    return output[:,:]
end
