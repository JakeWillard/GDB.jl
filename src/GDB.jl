module GDB

export test_simulation, Variable, Density, Temp
export jacobi, laplacian2d

using Distributed
@everywhere using LinearAlgebra
@everywhere include("./variables.jl")
@everywhere include("./jacobi.jl")
@everywhere include("./laplacian2d.jl")

# have placeholder test_simulation for now
function test_simulation()

    return true
end


end
