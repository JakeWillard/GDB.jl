module GDB

export test_simulation

using Distributed
using SharedArrays
using LinearAlgebra
using SparseArrays
using HDF5
using ForwardDiff
using RecipesBase

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




# have placeholder test_simulation for now
function test_simulation()

    return true
end


end
