module GDB

export test_simulation, Variable, Density, Temp

using Distributed
@everywhere using LinearAlgebra
@everywhere include("./variables.jl")

# have placeholder test_simulation for now
function test_simulation()

    return true
end


end
