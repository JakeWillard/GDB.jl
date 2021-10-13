

module DistributedSolver

export dlinsolve, prbgs

using Distributed
@everywhere using LinearAlgebra, SparseArrays

include("./pgmres.jl")
include("./prbgs.jl")

function dlinsolve(A::SparseMatrixCSC, x0::Vector{Float64}, b::Vector{Float64})


end


end
