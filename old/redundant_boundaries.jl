
function reflection_matrix(delta, mx, my, MinvT, bar::Barrier, grd::Grid)

    # define function we want to map onto grd
    function f(x, y, grd)
        x, y = reflection(x, y, delta, bar)
        row_dat, row_js = interpolation_row(x, y, mx, my, MinvT, grd)
        return hcat(row_dat, row_js)
    end

    M = grid_map(f, (mx*my,2), grd)
    dat = M[:,1]
    js = Int64[M[:,2]...]
    is = vcat([k*ones(Int64, mx*my) for k=1:grd.Nk]...)

    R2d = sparse(is, js, dat, grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
    R3d = kron(sparse(I, grd.Nz, grd.Nz), R2d)

    return R3d
end


function boundary_operators(mx, my, MinvT, deltas::Vector{Float64}, bars::Vector{Barrier}, qs::Vector{Float64}, grd::Grid)

    Nb = length(bars)

    PEN = SparseMatrixCSC[sparse(I, grd.Nz*grd.Nk, grd.Nz*grd.Nk)]
    REF = SparseMatrixCSC[]
    DCHLT = SparseMatrixCSC[]
    NMANN = SparseMatrixCSC[]

    for i=1:Nb
        P = Diagonal(f_to_grid((x,y)-> smoothstep(x, y, deltas[i], bars[i]), grd))
        for j=1:length(PEN)
            PEN[j] = P * PEN[j]
        end
        append!(PEN, [(I - P)/qs[i]])
    end

    for i=1:Nb
        R = reflection_matrix(deltas[i], mx, my, MinvT, bars[i], grd)
        @info "Calculated reflection matrix $i"
        append!(REF, [R])
        append!(DCHLT, [qs[i]^2 * 0.5*PEN[i+1]*(I + R)])
        append!(NMANN, [0.5*PEN[i+1]*(I - R)])
    end

    return PEN, REF, DCHLT, NMANN
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
