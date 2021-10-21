
module ImmersedMirrors


export Grid, coarse_grid, function_to_grid, operator_to_grid
export Mirror, distance_to_mirror, mirror_image, smoothstep
export Extrapolator, make_dirichlet

using Distributed, RecipesBase
@everywhere using LinearAlgebra, SparseArrays, ProgressMeter, ForwardDiff

include("./mirror.jl")
include("./grid.jl")
include("./extrapolator.jl")

end
