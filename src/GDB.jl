module GDB


using Distributed

# @everywhere begin
using DistributedArrays
using SharedArrays
using LinearAlgebra
using IterativeSolvers
using SparseArrays
using HDF5
using ForwardDiff
using RecipesBase
using Logging
using ProgressMeter

include("./domain/grid.jl")
include("./domain/barrier.jl")
include("./domain/interpolation.jl")
include("./domain/operators.jl")
include("./domain/ghosts.jl")

include("./linsolve/gmres.jl")
include("./linsolve/relaxation.jl")
include("./linsolve/solvers.jl")

include("./magflux/trace.jl")
include("./magflux/maps.jl")
include("./magflux/annalus.jl")
include("./magflux/solovev.jl")

include("./physics/formulas.jl")
include("./physics/timestep.jl")

# end



end
