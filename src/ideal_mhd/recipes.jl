

@userplot PlotSetup

@recipe function f(plt::PlotSetup)

    grd = plt.args[1].grd
    seriestype := :scatter
    markersize := 1

    @series begin
        grd.points[1,:], grd.points[2,:]
    end
end
