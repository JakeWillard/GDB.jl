

function checkerboard_function(x, y, h)

    i = Int(floor(x/h))
    j = Int(floor(y/h))

    a = mod(i, 2) == 0 ? -1 : 1
    b = mod(j, 2) == 0 ? -1 : 1

    return a*b
end


function checkerboard_setup(Nx, Ny)

    chain = []
    ys = LinRange(-1.5, 1.5, 300)
    xs = 0.75 .+ 0.1*sin.(ys*pi)
    chain = [[xs[i], ys[i]] for i=1:300]
    chain = [chain; [[-0.5, 2.0], [-0.5, -2.0]]]
    chain = hcat(chain...)
    M = Mirror([chain])

    grd = Grid([0.0, 0.0], [1.0, 1.0], Nx, Ny, 1) do x,y
        distance_to_mirror(x, y, M)[1] > -0.1
    end

    gd = GhostData(M, grd)

    fid = h5open("checkerboard_setup.h5", "w")
    save_grid(fid, grd, "Grid")
    save_ghost_data(fid, gd, "gd")
    fid["chain"] = chain[:,:]
    close(fid)
end


function checkerboard_plot(h)

    fid = h5open("checkerboard_setup.h5", "r")
    grd = load_grid(fid, "Grid")
    gd = load_ghost_data(fid, "gd")
    chain = fid["chain"][:,:]
    close(fid)

    fvec = f_to_grid(grd) do x,y
        sin(x*pi*10)
        # checkerboard_function(x, y, h)
    end

    p = heatmap(gd.R*fvec, grd, 1, aspect_ratio=:equal, color=:pastel, colorbar=false)
    plot!(p, chain[1,:], chain[2,:], axis=[], bordercolor="white", linewidth=3, color=:red, legend=false, xlims=(grd.r0[1], grd.r1[1]), ylims=(grd.r0[2], grd.r1[2]))
    return p
end


function plot_mirror_transform(N, h)

    fid = h5open("checkerboard_setup.h5", "r")
    grd = load_grid(fid, "Grid")
    gd = load_ghost_data(fid, "gd")
    chain = fid["chain"][:,:]
    close(fid)

    xs = zeros(2, grd.Nk)
    ys = zeros(2, grd.Nk)
    is, js, dats = findnz(gd.R)
    n = 0
    for k=1:grd.Nk
        x0, y0 = grd.points[:,is[k]]
        xf, yf = grd.points[:,js[k]]
        if norm([x0, y0] - [xf, yf]) > h
            n += 1
            xs[:,n] = [x0, xf]
            ys[:,n] = [y0, yf]
        end
    end

    skip = maximum([div(n, N), 1])
    xs = xs[:,1:skip:n]
    ys = ys[:,1:skip:n]

    p = plot(chain[1,:], chain[2,:], color=:red, linewidth=3)
    scatter!(p, grd.points[1,:], grd.points[2,:], color=:blue, markersize=1)
    plot!(p, xs, ys, color=:black, axis=[], bordercolor="white", arrow=true, linewidth=0.5, aspect_ratio=:equal, legend=false, xlims=(grd.r0[1], grd.r1[1]), ylims=(grd.r0[2], grd.r1[2]))
    return p
end


function old_code_grid(Nr, Nt)

    ts = LinRange(0, 2*pi, Nt+1)
    rs = LinRange(0.75, 1, Nr)
    points = zeros(2, Nr*Nt)
    k = 0
    for i=1:Nr
        for j=1:Nt
            k += 1
            points[:,k] = rs[i]*[cos(ts[j]), sin(ts[j])]
        end
    end
    inner = hcat(points[:,1:Nt], points[:,1])
    outer = hcat(points[:,end-Nt+1:end],points[:,end-Nt+1])

    p = plot(inner[1,:], inner[2,:])
    plot!(p, outer[1,:], outer[2,:])
    scatter!(points[1,:], points[2,:], axis=[], bordercolor="white", legend=false, aspect_ratio=:equal, markersize=2, color=:black)
    return p
end
