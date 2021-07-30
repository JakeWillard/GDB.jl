

# functions for derivative matrices
function stencil1d(m::Int)

    ic = Int(ceil(m/2.0))
    M = zeros((m, m))
    for i=1:m
        for j=1:m
            M[i, j] = (i - ic)^(j-1) / factorial(j-1)
        end
    end
    return inv(M)
end


# NOTE: Minv is the output of running stencil1d()
function regular_derivative_1d(N, n, Minv)

    D = sparse(zeros(Float64, (N, N)))
    stencil = Minv[n+1,:]
    m = length(stencil)
    ic = Int(ceil(m/2.0))

    for i=1:m
        k = i - ic
        vec = stencil[i] * ones(N - abs(k))
        D += spdiagm(k => vec)
    end

    return D
end


function derivative_matrix(nx, ny, Mxinv, Myinv, grd::Grid)

    Dx = regular_derivative_1d(grd._Nx, nx, Mxinv)
    Dy = regular_derivative_1d(grd._Ny, ny, Myinv)

    D = kron(Dy, Dx)
    P = grd.Proj
    Pinv = transpose(grd.Proj)

    return P * D * Pinv / ((grd.dx)^nx *(grd.dy)^ny)
end


function reflection_matrix(delta, mx, my, MinvT, bar::Barrier, grd::Grid)

    # define function we want to map onto grd
    function f(x, y, grd)
        x, y = reflection(x, y, delta, bar)
        row_dat, row_js = interpolation_row(x, y, mx, my, MinvT, grd)
        return hcat(row_dat, row_js, dS_vec)
    end

    M = grid_map(f, (grd.mx*grd.my,2), grd)
    dat = M[:,1]
    js = Int64[M[:,2]...]
    is = vcat([k*ones(Int64, grd.mx*grd.my) for k=1:Nk]...)

    return sparse(is, js, dat, Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
end


function boundary_operators(deltas::Vector{Float64}, bars::Vector{Barrier}, qs::Vector{Float64}, grd::Grid)

    Nk = size(grd.points)[2]
    Nb = length(bars)

    PEN = SparseMatrixCSC[sparse(I, Nk, Nk)]
    REF = SparseMatrixCSC[]
    DCHLT = SparseMatrixCSC[]
    NMANN = SparseMatrixCSC[]

    for i=1:Nb
        P = Diagonal(f_to_vec((x,y)-> smoothstep(x, y, delta[i], bars[i]), grd))
        PEN = P .* ops
        append!(ops, [(I - P)/qs[i]])
    end

    for i=1:Nb
        R = reflection_matrix(deltas[i], mx, my, MinvT, bars[i], grd)
        append!(REF, [R])
        append!(DCHLT, [0.5*(I + R)])
        append!(NMANN, [0.5*(I - R)])
    end

    return PEN, REF, PEN[2:end] .* DCHLT, PEN[2:end] .* NMANN
end
