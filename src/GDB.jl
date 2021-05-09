module GDB

export test_simulation

using Distributed
using LinearAlgebra
using SparseArrays

include("./prim/splines.jl")
include("./prim/grids.jl")
include("./prim/matrices.jl")
include("./prim/jacobi.jl")
include("./prim/data.jl")



# have placeholder test_simulation for now
function test_simulation()

    return true
end


end
