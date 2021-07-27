

# NOTE: tracing parallel or anti-parallel to b can be achieved with positive or negative ds
function fieldline_map_matrix(bx::Function, by::Function, bz::Function, ds::Float64, mx::Int64, my::Int64, MinvT::Matrix{Float64}, Nz::Int64, grd::Grid)

    deltaPhi = 2*pi / Nz

    # define function we want to map onto grd
    function f(x0, y0, grd)

        x, y, deltaS = trace_fieldline(x0, y0, bx, by, bz, ds, deltaPhi)
        row_dat, row_js = interpolation_row(x, y, mx, my, MinvT, grd)
        dS_vec = zeros(Float64, mx*my)
        dS_vec[1] = deltaS

        return hcat(row_dat, row_js, dS_vec)
    end

    M = grid_map(f, (grd.mx*grd.my,3), grd)
    dS = M[1:grd.mx*grd.my:end,3]
    dat = M[:,1]
    js = Int64[M[:,2]...]
    is = vcat([k*ones(Int64, grd.mx*grd.my) for k=1:grd.Nk]...)

    matrix = sparse(is, js, dat, grd.Nk, grd.Nx*grd.Ny) * transpose(grd.Proj)
    return matrix, dS
end


function fieldline_derivatives(bx::Function, by::Function, bz::Function, ds::Float64, mx::Int64, my::Int64, MinvT::Matrix{Float64}, Nz::Int64, grd::Grid)

    # compute forward and backward maps
    FM, dS1 = fieldline_map_matrix(bx, by, bz, ds, mx, my, MinvT, Nz, grd)
    BM, dS2 = fieldline_map_matrix(bx, by, bz, -ds, mx, my, MinvT, Nz, grd)

    # make diagonal matrices for distance coefficients
    _A = Diagonal(-dS1 ./ (dS2.^2 + dS1.*dS2))
    _B = Diagonal(1 ./ dS2 - 1 ./ dS1)
    _C = Diagonal(dS2 ./ (dS1.^2 + dS1.*dS2))
    _D = Diagonal(1 ./ (dS2.^2 + dS1.*dS2))
    _E = Diagonal(-1 ./ (dS1.*dS2))
    _F = Diagonal(1 ./ (dS1.^2 + dS1.*dS2))

    # make off diagonal matrices for kronecker product
    Ia = spdiag(-1 => ones(Float64, Nz-1))
    Ib = sparse(I, Nz, Nz)
    Ic = spdiag(1 => ones(Float64, Nz-1))

    # construct derivative matrices
    Ds = kron(Ia, _A*BM) + kron(Ib, _B) + kron(Ic, _C*FM)
    Dss = kron(Ia, _D*BM) + kron(Ib, _E) + kron(Ic, _F*FM)

    return Ds, Dss
end