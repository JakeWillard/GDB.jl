module GDB

export test_simulation

using Distributed
using LinearAlgebra
using SparseArrays


include("./prim/grids.jl")
include("./prim/penalize.jl")
include("./prim/jacobi.jl")
include("./prim/physical.jl")


# have placeholder test_simulation for now
function test_simulation()

    return true
end


end
