
module ImmersedMirrors


export Grid, function_to_grid, operator_to_grid
export Mirror, distance_to_mirror, mirror_image, smoothstep
export GhostData, flip_segments, constrain_system

using Distributed, RecipesBase
@everywhere using LinearAlgebra, SparseArrays, ProgressMeter, ForwardDiff

include("./mirrors.jl")
include("./grid.jl")
include("./ghosts.jl")

end
