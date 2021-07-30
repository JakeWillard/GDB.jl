using Distributed
using ForwardDiff
using SharedArrays
using LinearAlgebra
using SparseArrays
using HDF5
using RecipesBase
using Plots

include("/home/jake/Documents/GDB.jl/src/domain/grid.jl")
include("/home/jake/Documents/GDB.jl/src/domain/barrier.jl")
include("/home/jake/Documents/GDB.jl/src/domain/interpolation.jl")
include("/home/jake/Documents/GDB.jl/src/domain/operators.jl")

include("/home/jake/Documents/GDB.jl/src/linsolve/gmres.jl")
include("/home/jake/Documents/GDB.jl/src/linsolve/relaxation.jl")
include("/home/jake/Documents/GDB.jl/src/linsolve/solvers.jl")

include("/home/jake/Documents/GDB.jl/src/magflux/trace.jl")
include("/home/jake/Documents/GDB.jl/src/magflux/maps.jl")
include("/home/jake/Documents/GDB.jl/src/magflux/annalus.jl")

include("/home/jake/Documents/GDB.jl/src/physics/formulas.jl")
include("/home/jake/Documents/GDB.jl/experimental/calcs.jl")
include("/home/jake/Documents/GDB.jl/experimental/wrkgen.jl")



n0 = 1e13
T0 = 20
beta0 = 0.3
B0 = 22e3
q = 3
Nx = 50
Nz = 5

a = 67
R0 = 167
k = 0.2
dt = 0.001

wrk = rectangle_workspace(Nx, Nz, q, a, R0, k, n0, T0, B0, dt)

include("/home/jake/Documents/GDB.jl/experimental/recipes.jl")
