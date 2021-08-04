

function cluster_points(centers::Matrix{Float64}, grd::Grid)

    Nc = size(centers)[2]
    M = zeros(2, grd.Nk, Nc)
    J = zeros(Int64, (grd.Nk, Nc))
    ks = zeros(Int64, Nc)

    for i=1:grd.Nk
        j = sortperm([norm(grd.points[:,i] - centers[:,k]) for k=1:Nc])[1]
        ks[j] += 1
        M[:,ks[j],j] = grd.points[:,i]
        J[ks[j],j] = i
    end

    points = zeros(2, grd.Nk)
    is = Int64[1:grd.Nk...]
    js = zeros(Int64, grd.Nk)
    dat = ones(Int64, grd.Nk)
    k = 0
    for i=1:Nc
        for j=1:ks[i]
            k += 1
            points[:,k] = M[:,j,i]
            js[k] = J[j,i]
        end
    end

    return points, sparse(is, js, dat, grd.Nk, grd.Nk)
end


function sort_by_angle(grd::Grid, r0)

    thetas = Float64[atan(grd.points[2,i]-r0[2], grd.points[1,i]-r0[1]) for i=1:grd.Nk]

    is = Int64[1:grd.Nk...]
    js = sortperm(thetas)
    dat = ones(Int64, grd.Nk)

    return grd.points[:,js], sparse(is, js, dat, grd.Nk, grd.Nk)
end
