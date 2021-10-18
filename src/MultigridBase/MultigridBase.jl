
module MultigridBase

export AffineMap, LinearProblem, compute_residual!, coarsen_problem!, jsmooth!
export cartesian_1d_crs, cartesian_2d_crs

using LinearAlgebra, SparseArrays

include("./base.jl")
include("./cartesian.jl")


end
