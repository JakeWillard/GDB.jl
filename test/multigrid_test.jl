using Plots
using RecipesBase
using LinearAlgebra
using SparseArrays
using Distributed

include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/prim/walls.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/prim/grids.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/prim/matrices.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/prim/splines.jl")
include("C:/Users/lucas/OneDrive/Documents/Dartmouth/Research/Code/GDB.jl/src/prim/projection.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/prim/multigrid.jl")
include("C:/Users/lucas/OneDrive/Documents/GitHub/GDB.jl/src/ideal_mhd/setup.jl")




# derivatives from matrices without projections
# call simple to put into multigrid_grids
# have v_cycle_general spit out residual after one time through
# norm of residual divided by norm of rhs
# return x and E_fine, if norm(E_fine)/norm(b) stop
# A - laplacian
# b - f_to_grid, first argument is function (say, (x,y) -> sin(x) ), second argument is grid
# vec to mesh spits out vector that you can apply a heatmap to using above notation
# if singular, b = transpose(A) * b, A = transpose(A) * A

grid, wall, deltas, corners = simple(60, 60, 10, 1, 20, 20, 3)
multiresgrids = multigrid_grids(grid, [wall], deltas, corners, [3,3])
A = laplacian(grid)
b = f_to_grid((x,y) -> sin(x), grid)
x = A \ b


plotvar(x, grid, :heatmap)
