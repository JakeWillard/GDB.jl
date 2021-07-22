module GDB

export test_simulation

using Distributed
using SharedArrays
using LinearAlgebra
using SparseArrays
using HDF5
using RecipesBase

include("./prim/trace.jl")
include("./prim/walls.jl")
include("./prim/splines.jl")
include("./prim/grids.jl")
include("./prim/matrices.jl")
include("./prim/relaxation.jl")
include("./prim/gmres.jl")
include("./prim/simtools.jl")

include("./ideal_mhd/setup.jl")
include("./ideal_mhd/io.jl")
include("./ideal_mhd/recipes.jl")
include("./ideal_mhd/integration.jl")
include("./ideal_mhd/simulation.jl")





# have placeholder test_simulation for now
function test_simulation()

    return true
end


end
