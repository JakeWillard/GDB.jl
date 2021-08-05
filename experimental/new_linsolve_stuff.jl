
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




# struct BoundaryCondition
#
#     Proj :: SparseMatrixCSC
#     R :: SparseMatrixCSC
#
# end


struct GhostConditions

    Proj :: SparseMatrixCSC
    R :: SparseMatrixCSC
    ghost_indices

end


function GhostConditions(mx, my, MinvT, parities::Vector{Int64}, bars::Vector{Barrier}, grd::Grid)

    is = vcat([k*ones(Int64, mx*my) for k=1:grd.Nk]...)

    js = Int64[]
    dat = Float64[]
    j_proj = Int64[]

    ghost_indices = zeros(Int64, grd.Nk, length(bars))
    ng = zeros(Int64, length(bars))

    for k=1:grd.Nk
        us = [smoothstep(grd.points[:,k]..., 1.0, bar) for bar in bars]
        if !isempty(us[us .< 0.5])

            us[us .> 0.5] .= -Inf
            bind = sortperm(us)[end]
            barr = bars[bind]
            par = parities[bind]

            ng[bind] += 1
            ghost_indices[ng[bind], bind] = k

            x, y = barr.rmap(grd.points[:,k]...)
            row_dat, row_j = interpolation_row(x, y, mx, my, MinvT, grd)
            js = [js; row_j]
            dat = [dat; par*row_dat]
            j_proj = [j_proj; k]
        else

            row_dat, row_j = interpolation_row(grd.points[:,k]..., mx, my, MinvT, grd)
            js = [js; row_j]
            dat = [dat; row_dat]
        end
    end

    Proj = sparse(Int64[1:length(j_proj)...], j_proj, ones(length(j_proj)), length(j_proj), grd.Nk)
    R = sparse(is, js, dat, grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)

    return GhostConditions(Proj, R, ghost_indices)
end


function swap_parity(gc::GhostConditions, inds...)

    S = I + zeros(size(gc.R))

    for i in inds
        gi = gc.ghost_indices[:,i]
        gi = gi[gi .> 0]
        S[gi,:] = -S[gi,:]
    end

    R = S * gc.R
    return GhostConditions(gc.Proj, R, gc.ghost_indices)
end





# function BoundaryCondition(condition_type::Int64, mx, my, MinvT, orientation::Int64, bar::Barrier, grd::Grid)
#
#     is_r = Int64[]
#     js_r = Int64[]
#     dat_r = Float64[]
#
#     is_p = Int64[1:grd.Nk...]
#     js_p = Int64[]
#     dat_p = ones(Float64, grd.Nk)
#
#     if condition_type == 1
#         parity = -1   # Dirichlet-like condition
#     elseif condition_type == 2
#         parity = 1    # Neumann-like condition
#     else
#         @error "Invalid condition type"
#     end
#
#     ng = 0
#     for k=1:grd.Nk
#         if smoothstep(grd.points[:,k]..., orientation, bar) < 0.5
#             ng += 1
#             x, y = bar.rmap(grd.points[:,k]...)
#             row_dat, row_j = interpolation_row(x, y, mx, my, MinvT, grd)
#
#             is_r = [is_r; [k for _=1:mx*my]]
#             js_r = [js_r; row_j]
#             dat_r = [dat_r; parity * row_dat]
#
#             js_p = [js_p; k]
#         else
#             row_dat, row_j = interpolation_row(grd.points[:,k]..., mx, my, MinvT, grd)
#             is_r = [is_r; [k for _=1:mx*my]]
#             js_r = [js_r; row_j]
#             dat_r = [dat_r; row_dat]
#         end
#
#     end
#     @info "Number of ghosts:" ng
#
#     R = sparse(is_r, js_r, dat_r, grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
#     P = sparse(is_p[1:ng], js_p, dat_p[1:ng], ng, grd.Nk)
#
#     return BoundaryCondition(P, R)
# end
