

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

    grd = Grid([1-(1/1.618), 0.0], [1.0, 1.0], Nx, Ny, 1) do x,y
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

    p = heatmap(gd.R*fvec, grd, 1, aspect_ratio=:equal, color=:greys, colorbar=false)
    plot!(p, chain[1,:], chain[2,:], axis=[], bordercolor="white", linewidth=3, color=:red, size=(1000, 1618), legend=false, xlims=(grd.r0[1], grd.r1[1]), ylims=(grd.r0[2], grd.r1[2]))
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

    p = plot(chain[1,:], chain[2,:], color=:red, linewidth=3, )
    scatter!(p, grd.points[1,:], grd.points[2,:], color=:blue, markersize=3, )
    plot!(p, xs, ys, color=:black, axis=[], markersize=2, bordercolor="white", arrow=true, size=(1000, 1618), linewidth=0.5, aspect_ratio=:equal, legend=false, xlims=(grd.r0[1], grd.r1[1]), ylims=(grd.r0[2], grd.r1[2]))
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


function new_code_grid(Nx, Ny)

    boundary_setup("./poster_boundary_setup.h5", 40, 40, 40)

    fid = h5open("./poster_boundary_setup.h5", "r")
    M = load_mirror(fid, "Mirror")
    inner_boundary = fid["inner_boundary"][:]
    outer_boundary = fid["outer_boundary"][:]
    close(fid)

    grd = Grid(Nx, Ny, 1) do
        inside(x,y) = distance_to_mirror(x, y, M)[1] > -0.05
        r0 = Float64[0.5, -1]
        r1 = Float64[1.5, 1]
        inside, r0, r1
    end

    c = 0
    inner = []
    while true
        c += 1
        append!(inner, [M.verts[:,c]])
        if norm(M.verts[:,c] - M.verts[:,c+1]) > 0.2
            break
        end
    end
    inner = hcat(inner..., inner[1])
    outer = hcat(M.verts[:,c+1:end], M.verts[:,c+1])

    p = plot(inner[1,:], inner[2,:])
    plot!(p, outer[1,:], outer[2,:])
    scatter!(p, grd.points[1,:], grd.points[2,:], axis=[], bordercolor="white", legend=false, aspect_ratio=:equal, markersize=2, color=:black)
    return p
end


function setup_advection_diffusion(Nx, Ny)

    boundary_setup("./poster_boundary_setup.h5", 40, 40, 40)

    fid = h5open("./poster_boundary_setup.h5", "r")
    M = load_mirror(fid, "Mirror")
    inner_boundary = fid["inner_boundary"][:]
    outer_boundary = fid["outer_boundary"][:]
    close(fid)

    grd = Grid(Nx, Ny, 1) do
        inside(x,y) = distance_to_mirror(x, y, M)[1] > -0.05
        r0 = Float64[0.5, -1]
        r1 = Float64[1.5, 1]
        inside, r0, r1
    end

    gd = GhostData(M, grd)
