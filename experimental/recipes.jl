

@userplot PlotVar
@recipe function f(pv::PlotVar)

    var, index, grd, st = pv.args
    z = vec_to_mesh(var, grd)[:,:,index]
    seriestype := st

    @series begin
        transpose(z)
    end
end


@userplot PlotGrid
@recipe function f(pg::PlotGrid)

    grd, clr = pg.args
    z = vec_to_mesh(ones(Float64, grd.Nk*grd.Nz), grd)[:,:,1]
    seriestype := :heatmap
    fillcolor := clr
    colorbar := false
    @series begin transpose(z) end
end


# @userplot PlotSetup
# @recipe function f(ps::PlotSetup)
#
#     stp = ps.args[1]
#     x = stp.grd.points[1,:]
#     y = stp.grd.points[2,:]
#     pvec = diag(stp.P1*stp.P2*stp.P3)
#
#     p0 = vec_to_mesh(pvec, stp.grd)
#     p1 = vec_to_mesh(stp.R1*pvec, stp.grd)
#     p2 = vec_to_mesh(stp.R2*pvec, stp.grd)
#     p3 = vec_to_mesh(stp.R3*pvec, stp.grd)
#
#     layout := @layout [b c; d e]
#     seriestype := :heatmap
#     colorbar := false
#
#     @series begin
#         subplot := 1
#         p0
#     end
#
#     @series begin
#         subplot := 2
#         p1
#     end
#
#     @series begin
#         subplot := 3
#         p2
#     end
#
#     @series begin
#         subplot := 4
#         p3
#     end
# end
