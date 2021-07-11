

@userplot PlotVar
@recipe function f(pv::PlotVar)

    var, grd, st = pv.args
    z = vec_to_mesh(var, grd)
    seriestype := st

    @series begin
        z
    end
end


@userplot PlotGrid
@recipe function f(pg::PlotGrid)

    grd = pg.args[1]
    x = grd.points[1,:]
    y = grd.points[2,:]
    seriestype := :scatter
    markersize := 1

    @series begin
        x, y
    end
end
