module GDB

export test_simulation, Variable, Density, Temp
export jacobi, laplacian2d, dist_jacobi

using Distributed
@everywhere using LinearAlgebra
@everywhere using SparseArrays
@everywhere include("./variables.jl")
@everywhere include("./jacobi.jl")
@everywhere include("./laplacian2d.jl")
@everywhere include("./p_jacobi.jl")

# have placeholder test_simulation for now
function test_simulation()

    return true
end


end
