

@userplot PlotVar
@recipe function f(pv::PlotVar)

    var, grd, st = pv.args
    z = vec_to_mesh(var, grd)
    seriestype := st

    @series begin
        z
    end
end


@userplot PlotSetup
@recipe function f(ps::PlotSetup)

    stp = ps.args[1]
    x = stp.grd.points[1,:]
    y = stp.grd.points[2,:]
    pvec = diag(stp.P1*stp.P2*stp.P3)

    # p0 = zeros(Float64, stp.grd.Nx*stp.grd.Ny)
    # p1 = zeros(Float64, stp.grd.Nx*stp.grd.Ny)
    # p2 = zeros(Float64, stp.grd.Nx*stp.grd.Ny)
    # p3 = zeros(Float64, stp.grd.Nx*stp.grd.Ny)
    #
    # p0[:] = pvec
    # p1[:] = stp.R1*pvec
    # p2[:] = stp.R2*pvec
    # p3[:] = spt.R3*pvec
    #
    # p0_mesh = vec_to_mesh(stp.R1*)

    p0 = vec_to_mesh(pvec, stp.grd)
    p1 = vec_to_mesh(stp.R1*pvec, stp.grd)
    p2 = vec_to_mesh(stp.R2*pvec, stp.grd)
    p3 = vec_to_mesh(stp.R3*pvec, stp.grd)

    layout := @layout [a _; b c; d e]

    @series begin
        subplot := 1
        seriestype := :scatter
        markersize := 1
        x, y
    end

    seriestype := :heatmap
    colorbar := false

    @series begin
        subplot := 2
        p0
    end

    @series begin
        subplot := 3
        p1
    end

    @series begin
        subplot := 4
        p2
    end

    @series begin
        subplot := 5
        p3
    end
end
