using GDB
using Test
using SparseArrays
using LinearAlgebra
using Plots


function test_domain()

    bar1 = GDB.Barrier() do
        func(x,y) = x
        val = 1.0
        rmap(x,y) = 1 - x, y
        func, val, rmap, -1
    end

    bar2 = GDB.PolyBarrier() do
        Float64[1 1; -1 2]
    end

    bar3 = GDB.Everywhere()
    bar4 = GDB.Nowhere()

    grd = GDB.Grid(50, 50, 1) do
        r0 = Float64[0, 0]
        r1 = Float64[2, 1]
        inside(x,y) = GDB.check_if_inside(x, y, [0.1], [bar1])
        inside, r0, r1
    end

    fvec = GDB.f_to_grid(grd) do x,y
        sin(x)*sin(y)
    end

    h = heatmap((fvec, grd, 1))

    f_func = GDB.interpolation_function(fvec, 2, 2, GDB.stencil2d(2, 2), grd)

    gc = GDB.GhostConditions(2, 2, GDB.stencil2d(2, 2), [bar1], grd)
    gc2 = GDB.swap_sign(gc, 1)

    A = sparse(I, grd.Nk, grd.Nk)
    xb = zeros(grd.Nk)
    x1 = rand(grd.Nk)
    x = gc.Proj * x1

    x2 = GDB.extrapolate_ghosts(x, xb, gc)
    Anew, B = GDB.require_boundary_conditions(A, gc)

    Minv = GDB.stencil1d(3)
    Dx = GDB.derivative_matrix(1, 0, Minv, Minv, grd)

    avg = GDB.line_average(Float64[0.5 1.5; 0.5 0.5], 2, 2, GDB.stencil2d(2, 2), grd)

    return true
end




@testset "domain" begin
    @test test_domain()
end
