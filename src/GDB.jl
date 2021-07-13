module GDB

export test_simulation

using Distributed
using LinearAlgebra
using SparseArrays
using HDF5
using RecipesBase

include("/home/jake/Documents/GDB.jl/src/prim/trace.jl")
include("/home/jake/Documents/GDB.jl/src/prim/walls.jl")
include("/home/jake/Documents/GDB.jl/src/prim/splines.jl")
include("/home/jake/Documents/GDB.jl/src/prim/grids.jl")
include("/home/jake/Documents/GDB.jl/src/prim/matrices.jl")
include("/home/jake/Documents/GDB.jl/src/prim/vars.jl")

include("/home/jake/Documents/GDB.jl/src/ideal_mhd/setup.jl")
include("/home/jake/Documents/GDB.jl/src/ideal_mhd/io.jl")
include("/home/jake/Documents/GDB.jl/src/ideal_mhd/recipes.jl")



# have placeholder test_simulation for now
function test_simulation()

    return true
end


end
