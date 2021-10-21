

mutable struct MultigridSolver{T}

    p :: Vector{LinearProblem{T}}  # linear problems
    r :: Vector{AffineMap{T}}  # restrictors
    i :: Vector{AffineMap{T}}   # interpolators
    n :: Int64

end
function MultigridSolver{T}(L::LinearProblem{T}, max_layers::Int64) where {T}
    p = Vector{LinearProblem{T}}(undef, max_layers)
    r = Vector{AffineMap{T}}(undef, max_layers-1)
    i = Vector{AffineMap{T}}(undef, max_layers-1)

    p[1] = L
    return MultigridSolver{T}(p, r, i, 1)
end
MultigridSolver(args...) = MultigridSolver{Float64}(args...)


function add_layer!(M::MultigridSolver{T}, Restr::AffineMap{T}, Interp::AffineMap{T}) where {T}

    Lc = coarsen_problem!(M.p[M.n], Restr, Interp)
    M.p[M.n+1] = Lc
    M.r[M.n] = Restr
    M.i[M.n] = Interp
    M.n += 1
end


function vcycle!(M::MultigridSolver{T}, nj::Vector{Int64}) where {T}

    # initial smoothing
    jsmooth!(M.p[1], nj[1])

    for i=2:M.n
        M.p[i].x[:] = M.r[i-1] * M.p[i-1].x
        compute_residual!(M.p[i])
        jsmooth!(M.p[i], nj[i])
    end

    for i=M.n-1:1
        M.p[i].x[:] = M.i[i] * M.p[i+1].x
        compute_residual!(M.p[i])
        jsmooth!(M.p[i], nj[i])
    end
end
