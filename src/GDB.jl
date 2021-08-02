module GDB

export test_simulation

using Distributed
using DistributedArrays
using SharedArrays
using LinearAlgebra
using SparseArrays
using HDF5
using ForwardDiff
using RecipesBase
using Logging

include("./domain/grid.jl")
include("./domain/barrier.jl")
include("./domain/interpolation.jl")
include("./domain/operators.jl")

include("./linsolve/gmres.jl")
include("./linsolve/relaxation.jl")
include("./linsolve/solvers.jl")

include("./magflux/trace.jl")
include("./magflux/maps.jl")
include("./magflux/annalus.jl")
include("./magflux/solovev.jl")

include("./physics/formulas.jl")
include("./physics/parameters.jl")
include("./physics/timestep.jl")

include("./preprocess/io.jl")
include("./preprocess/rectangle.jl")


# have placeholder test_simulation for now
function test_simulation()

    return true
end


end
