

mutable struct MultigridSolver{T}

    problems :: Vector{LinearProblem{T}}
    restrictors :: Vector{AffineMap{T}}
    interpolators :: Vector{AffineMap{T}}
    _m :: Int64

end
function MultigridSolver{T}(L::LinearProblem{T}, max_layers::Int64) where {T}
    problems = Vector{LinearProblem{T}}(undef, max_layers)
    restrictors = Vector{AffineMap{T}}(undef, max_layers)
    interpolators = Vector{AffineMap{T}}(undef, max_layers)

    problems[1] = L
    return MultigridSolver{T}(problems, restrictors, interpolators, 2)
end
MultigridSolver(args...) = MultigridSolver{Float64}(args...)


function add_layer!(M::MultigridSolver{T}, Restr::AffineMap{T}, Interp::AffineMap{T}) where {T}

    Lc = coarsen_problem!(M.problems[_m-1], Restr, Interp)
    M.problems[_m] = Lc
    M.restrictors[_m] = Restr
    M.interpolators[_m] = Interp
    M._m += 1
end
