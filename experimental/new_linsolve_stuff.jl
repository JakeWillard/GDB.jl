
"""
changes:
    - reflections should be calculated exclusively on ghost points.
    - reflection_matrix() function should have input argument determining if we want
      a reflection or anti-reflection operation, and also input for what to do on other points.
"""


struct LinearMap

    M :: SparseMatrixCSC
    f :: Vector{Float64}

    Base.size(l::LinearMap) = length(l.f)
    Base.getindex(l::LinearMap, inds) = LinearMap(l.M[inds,:], l.f[inds])
    Base.:*(l::LinearMap, v::Vector{Float64}) = l.M * v + f
    LinearMap(A::SparseMatrixCSC, b::Vector{Float64}) = new(A, b)
end





struct GhostConditions

    Proj :: SparseMatrixCSC
    Extr :: SparseMatrixCSC
    swaps :: Matrix{Float64}

end


function GhostConditions(mx, my, MinvT, bars::Vector{Barrier}, grd::Grid)

    ks = [Int64[] for _=1:length(bars)+1]
    js = [Int64[] for _=1:length(bars)+1]
    dats = [Float64[] for _=1:length(bars)+1]
    swaps = ones(Float64, (grd.Nk, length(bars)))
    ng = 0

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

            # enforce a minimum displacement, otherwise we get intolerable artifacts at the grid scale.
            dxmin = sqrt(grd.dx^2 + grd.dy^2)
            dx = Float64[x, y] - grd.points[:,k]
            if norm(dx) < dxmin
                x += dxmin*dx[1]/norm(dx)
                y += dxmin*dx[2]/norm(dx)
            end

            # add reflection rows
            row_dat, row_j = interpolation_row(x, y, mx, my, MinvT, grd)
            js[l] = [js[l]; row_j]
            dats[l] = [dats[l]; row_dat]

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
    M = sparse(vcat([k*ones(2*mx*my) for k=1:ng]...), vcat(js...), vcat(dats...), ng, grd._Nx*grd._Ny) * transpose(grd.Proj) * transpose(Perm)

    # invert to get extrapolation operation, delete small values
    N = inv(Matrix(M[1:ng, 1:ng]))
    N[abs.(N) .< 1e-8] .= 0
    R = -sparse(N) * M[:, ng+1:end]

    # make extrapolation matrix
    Extr_perm = sparse(Diagonal(ones(Float64, grd.Nk)))
    Extr_perm = sparse([1:grd.Nk]..., [1:grd.Nk]..., ones(Float64, grd.Nk), grd.Nk, grd.Nk)
    Extr_perm[1:ng, ng+1:end] = R
    Extr_perm[1:ng, 1:ng] .= 0
    Extr = transpose(Perm) * Extr_perm * Perm

    # compute projection
    Proj = sparse([1:grd.Nk-ng...], ks[end], ones(grd.Nk-ng), grd.Nk-ng, grd.Nk)

    return GhostConditions(Proj, Extr, swaps)
end









#
#
#
#
# function GhostConditions(mx, my, MinvT, parities::Vector{Int64}, bars::Vector{Barrier}, grd::Grid)
#
#     # is = vcat([k*ones(Int64, mx*my) for k=1:grd.Nk]...)
#     is = Int64[]
#     js = Int64[]
#     dat = Float64[]
#     j_proj = Int64[]
#
#     ghost_indices = zeros(Int64, grd.Nk, length(bars))
#     ng = zeros(Int64, length(bars))
#
#     count = 0
#     for k=1:grd.Nk
#         us = [smoothstep(grd.points[:,k]..., 1.0, bar) for bar in bars]
#         if !isempty(us[us .< 0.5])
#
#             us[us .> 0.5] .= -Inf
#             bind = sortperm(us)[end]
#             barr = bars[bind]
#             par = parities[bind]
#
#             count += 1
#             ng[bind] += 1
#             ghost_indices[ng[bind], bind] = k
#
#             x, y = barr.rmap(grd.points[:,k]...)
#             row_dat, row_j = interpolation_row(x, y, mx, my, MinvT, grd)
#             is = [is; count*ones(Int64, mx*my)]
#             js = [js; row_j]
#             dat = [dat; par*row_dat]
#
#             row_dat, row_j = interpolation_row(grd.points[:,k]..., mx, my, MinvT, grd)
#             is = [is; count*ones(Int64, mx*my)]
#             js = [js; row_j]
#             dat = [dat; row_dat]
#         else
#
#             # row_dat, row_j = interpolation_row(grd.points[:,k]..., mx, my, MinvT, grd)
#             # is = [is; k*ones(Int64, mx*my)]
#             # js = [js; row_j]
#             # dat = [dat; row_dat]
#             j_proj = [j_proj; k]
#         end
#     end
#
#     Proj = sparse(Int64[1:length(j_proj)...], j_proj, ones(length(j_proj)), length(j_proj), grd.Nk)
#     R = sparse(is, js, dat, grd.Nk - length(j_proj), grd._Nx*grd._Ny) * transpose(grd.Proj)
#     return R
#
#     return GhostConditions(Proj, R, ghost_indices)
# end


function swap_parity(gc::GhostConditions, inds...)

    Extr = gc.Extr

    for i in inds
        Extr = Diagonal(gc.swaps[:,i]) * Extr
    end

    return GhostConditions(gc.Proj, Extr, gc.swaps)
end


function extrapolate_ghosts(x::Vector{Float64}, xb::Vector{Float64}, gc::GhostConditions)

    return gc.Extr*transpose(gc.Proj)*x - (I - gc.Extr)*xb
end


function require_boundary_conditions(A::SparseMatrixCSC, gc::GhostConditions)

    P = gc.Proj
    Pt = transpose(P)

    Anew = P*A*gc.Extr*Pt
    Bterm = P*A*(I - gc.Extr)

    return Anew, Bterm
end


function require_boundary_conditions(L::LinearMap, xb::Vector{Float64}, gc::GhostConditions)

    Anew, Bterm = require_boundary_conditions(L.M, gc)
    fnew = gc.Proj*L.f + Bterm*xb

    return LinearMap(Anew, fnew)
end


function solve_pde(A::SparseMatrixCSC, b::Vector{Float64}, xb::Vector{Float64}, gc::GhostConditions)

    L = require_boundary_conditions(LinearMap(A, -b), xb, gc)
    x = L.M \ -L.f

    return extrapolate_ghosts(x, xb, gc)
end
