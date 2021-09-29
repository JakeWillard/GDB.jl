
module ImmersedMirrors


export Grid, function_to_grid, operator_to_grid
export Mirror, smoothstep
export GhostData, flip_segments, extrapolate_values, require_boundary_conditions

using Distributed
@everywhere using LinearAlgebra, SparseArrays, ProgressMeter, ForwardDiff
using Logging, RecipesBase

include("./mirrors.jl")
include("./grid.jl")
include("./ghosts.jl")

end
