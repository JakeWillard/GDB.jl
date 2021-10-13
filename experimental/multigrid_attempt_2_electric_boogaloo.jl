
function intergrid_1d_interp(N, m)

    # make stencil
    ic = Int(ceil(m/2.0))
    M = zeros((m, m))
    for i=1:m
        for j=1:m
            M[i, j] = (i - ic)^(j-1) / factorial(j-1)
        end
    end
    Minv = inv(M)
    sten = transpose(Minv) * Float64[(0.5)^(j-1)/factorial(j-1) for j=1:m]

    Ih = zeros(Float64, (N,N))
    Is = I + Ih
    for i=1:m
        k = i - ic
        vec = sten[i] * ones(N - abs(k))
        Ih += diagm(k => vec)
    end

    Iout = zeros(Float64, (2*N,N))
    Iout[1:2:end,:] = Is
    Iout[2:2:end,:] = Ih
    return Iout
end


struct TwoGrids

    coarsegrid :: Grid
    finegrid :: Grid
    Interp :: SparseMatrixCSC
    Restr :: SparseMatrixCSC

end


function TwoGrids(init::Function, p::Int64, m::Int64)

    coarsegrid = Grid(init()..., p-1; Nbuffer=0)
    finegrid = Grid(init()..., p; Nbuffer=0)

    compute cartesian restriction and interpolation matrices
    Ix = intergrid_1d_interp(coarsegrid._Nx, m)
    Iy = intergrid_1d_interp(coarsegrid._Ny, m)
    Ic = kron(sparse(Iy), sparse(Ix))
    Rc = transpose(Ic) / 4.0

    Interp = finegrid.Proj * Ic * transpose(coarsegrid.Proj)
    Restr = coarsegrid.Proj * Rc * transpose(finegrid.Proj)

    return TwoGrids(coarsegrid, finegrid, Interp, Restr)
end


function function_to_both_grids(f::Function, grds::TwoGrids)

    v_coarse = function_to_grid(f, grds.coarsegrid)
    v_fine = function_to_grid(f, grds.finegrid)
    v_coarse, v_fine
end


function operator_to_both_grids(O::SparseMatrixCSC, grds::TwoGrids)

    O_coarse = operator_to_grid(O, grds.coarsegrid)
    O_fine = operator_to_grid(O, grds.finegrid)
    O_coarse, O_fine
end


function both_extrapolators(M::Mirror, grds::TwoGrids)

    return Extrapolator(M, grds.coarsegrid), Extrapolator(M, grds.finegrid)
end





function two_level_vcycle(Ac, Af, Af_diag, Restr, Interp, x0, b, Nj)

    r = b - Af*x0
    Nf = length(x0)
    delta = zeros(Nf)

    # do fine smoothing
    for _=1:Nj
        delta[:] = r .\ Af_diag
        r[:] = r - Af*delta
    end

    # solve on coarse grid, interpolate back to coarse grid
    delta[:] = Interp * (Ac \ (Restr * r))
    r[:] = r - Af*delta

    # do one more fine smoothing
    for _=1:Nj
        delta = r .\ Af_diag
        r[:] = r - Af*delta
    end

    return x0 + delta
end
