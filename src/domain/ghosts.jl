


struct GhostConditions

    Proj :: SparseMatrixCSC
    Extr :: SparseMatrixCSC
    swaps :: Matrix{Float64}

end


function GhostConditions(mx, my, MinvT, bars::Vector{Barrier}, grd::Grid; w=0.5)

    ks = [Int64[] for _=1:length(bars)+1]
    js = [Int64[] for _=1:length(bars)+1]
    dats = [Float64[] for _=1:length(bars)+1]
    swaps = ones(Float64, (grd.Nk, length(bars)))
    ng = 0
    dr = 1.5*sqrt(grd.dx^2 + grd.dy^2)

    for k=1:grd.Nk
        us = [smoothstep(grd.points[:,k]..., 1.0, bar) for bar in bars]
        if !isempty(us[us .< 0.5])

            us[us .> 0.5] .= -Inf
            l = sortperm(us)[end]
            ng += 1

            # add permutation row
            ks[l] = [ks[l]; k]

            # compute reflection
            x, y = bars[l].rmap(grd.points[:,k]...)

            # for robustness, compute coefficients in terms of neighboring points.
            dist = Float64[x, y] - grd.points[:,k]
            er = dist / norm(dist)
            x1 = x + dr*er[1]
            x2 = x + 2*dr*er[1]
            y1 = y + dr*er[2]
            y2 = y + 2*dr*er[2]

            # add reflection rows
            row_dat1, row_j1 = interpolation_row(x1, y1, mx, my, MinvT, grd)
            row_dat2, row_j2 = interpolation_row(x2, y2, mx, my, MinvT, grd)
            js[l] = [js[l]; row_j1]
            js[l] = [js[l]; row_j2]
            dats[l] = [dats[l]; 2*row_dat1]
            dats[l] = [dats[l]; -1*row_dat2]

            # enforce a minimum displacement, otherwise we get intolerable artifacts at the grid scale.
            # dxmin = sqrt(grd.dx^2 + grd.dy^2)
            # dx = Float64[x, y] - grd.points[:,k]
            # if norm(dx) < dxmin
            #     x += dxmin*dx[1]/norm(dx)
            #     y += dxmin*dx[2]/norm(dx)
            # end

            # add reflection rows
            # row_dat, row_j = interpolation_row(x, y, mx, my, MinvT, grd)
            # js[l] = [js[l]; row_j]
            # dats[l] = [dats[l]; row_dat]

            # add identity rows
            row_dat, row_j = interpolation_row(grd.points[:,k]..., mx, my, MinvT, grd)
            js[l] = [js[l]; row_j]
            dats[l] = [dats[l]; row_dat]

            # change swaps
            swaps[k,l] = -1.0

        else

            # add permutation row
            ks[end] = [ks[end]; k]
        end
    end

    # compute permutation and reflection condition matrix
    Perm = sparse([1:grd.Nk...], vcat(ks...), ones(grd.Nk), grd.Nk, grd.Nk)
    M = sparse(vcat([k*ones(3*mx*my) for k=1:ng]...), vcat(js...), vcat(dats...), ng, grd._Nx*grd._Ny) * transpose(grd.Proj) * transpose(Perm)

    # invert to get extrapolation operation, delete small values
    N = inv(Matrix(M[1:ng, 1:ng]))
    N[abs.(N) .< 1e-8] .= 0
    R = -sparse(N) * M[:, ng+1:end]

    # make extrapolation matrix
    Extr_perm = sparse(Diagonal(ones(Float64, grd.Nk)))
    # Extr_perm = sparse([1:grd.Nk]..., [1:grd.Nk]..., ones(Float64, grd.Nk), grd.Nk, grd.Nk) # XXX
    Extr_perm[1:ng, ng+1:end] = R
    Extr_perm[1:ng, 1:ng] .= 0
    Extr = transpose(Perm) * Extr_perm * Perm

    # compute projection
    Proj = sparse([1:grd.Nk-ng...], ks[end], ones(grd.Nk-ng), grd.Nk-ng, grd.Nk)

    return GhostConditions(Proj, Extr, swaps)
end


function swap_sign(gc::GhostConditions, inds...)

    Extr = gc.Extr

    for i in inds
        Extr = Diagonal(gc.swaps[:,i]) * Extr
    end

    return GhostConditions(gc.Proj, Extr, gc.swaps)
end


function extrapolate_ghosts(x::Vector{Float64}, xb::Vector{Float64}, gc::GhostConditions)

    return gc.Extr*transpose(gc.Proj)*x + (I - gc.Extr)*xb
end


function require_boundary_conditions(A::SparseMatrixCSC, gc::GhostConditions)

    P = gc.Proj
    Pt = transpose(P)

    Anew = P*A*gc.Extr*Pt
    Bterm = P*A*(I - gc.Extr)

    return Anew, Bterm
end
