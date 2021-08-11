

function fieldline_images(bx::Function, by::Function, bz::Function, ds::Float64, grd::Grid)

    deltaPhi = 2*pi / grd.Nz

    img = DArray((6, grd.Nk), workers(), [1, length(workers())]) do inds

        Nk = length(inds[2])
        points0 = grd.points[:, inds[2]]
        imgpoints = zeros(6, Nk)

        for k=1:Nk

            imgpoints[1:3,k] = trace_fieldline(points0[:,k]..., bx, by, bz, ds, deltaPhi)
            imgpoints[4:6,k] = trace_fieldline(points0[:,k]..., bx, by, bz, -ds, deltaPhi)

        end
        imgpoints
    end

    return Matrix(img)
end


function fieldline_derivatives(img::Matrix{Float64}, mx::Int64, my::Int64, MinvT::Matrix{Float64}, grd::Grid)

    dS1 = img[3,:]
    dS2 = img[6,:]
    fwd_pts = img[1:2,:]
    bck_pts = img[4:5,:]

    # make diagonal matrices for distance coefficients
    _A = Diagonal(-dS1 ./ (dS2.^2 + dS1.*dS2))
    _B = Diagonal(1 ./ dS2 - 1 ./ dS1)
    _C = Diagonal(dS2 ./ (dS1.^2 + dS1.*dS2))
    _D = Diagonal(1 ./ (dS2.^2 + dS1.*dS2))
    _E = Diagonal(-1 ./ (dS1.*dS2))
    _F = Diagonal(1 ./ (dS1.^2 + dS1.*dS2))

    # make off diagonal matrices for kronecker product
    Ia = spdiagm(-1 => ones(Float64, grd.Nz-1))
    Ia[1,end] = 1
    Ib = sparse(I, grd.Nz, grd.Nz)
    Ic = spdiagm(1 => ones(Float64, grd.Nz-1))
    Ic[end,1] = 1

    # make fieldline map matrices
    is = vcat([k*ones(mx*my) for k=1:grd.Nk]...)
    fwd_j = Int64[]
    fwd_dat = Float64[]
    bck_j = Int64[]
    bck_dat = Float64[]
    for k=1:grd.Nk
        fwd_row_dat, fwd_row_j = interpolation_row(fwd_pts[:,k]..., mx, my, MinvT, grd)
        bck_row_dat, bck_row_j = interpolation_row(bck_pts[:,k]..., mx, my, MinvT, grd)
        fwd_j = [fwd_j; fwd_row_j]
        fwd_dat = [fwd_dat; fwd_row_dat]
        bck_j = [bck_j; bck_row_j]
        bck_dat = [bck_dat; bck_row_dat]
    end

    FWD = sparse(is, fwd_j, fwd_dat, grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
    BCK = sparse(is, bck_j, bck_dat, grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)

    # compute fieldline derivatives
    Ds = kron(Ia, _A*BCK) + kron(Ib, _B) + kron(Ic, _C*FWD)
    Dss = kron(Ia, _D*BCK) + kron(Ib, _E) + kron(Ic, _F*FWD)

    return Ds, Dss
end

#
#
# # NOTE: tracing parallel or anti-parallel to b can be achieved with positive or negative ds
# function fieldline_map_matrix(bx::Function, by::Function, bz::Function, ds::Float64, mx::Int64, my::Int64, MinvT::Matrix{Float64}, Nz::Int64, grd::Grid)
#
#     deltaPhi = 2*pi / Nz
#
#     # define function we want to map onto grd
#     function f(x0, y0, grd)
#
#         x, y, deltaS = trace_fieldline(x0, y0, bx, by, bz, ds, deltaPhi)
#         row_dat, row_js = interpolation_row(x, y, mx, my, MinvT, grd)
#         dS_vec = zeros(Float64, mx*my)
#         dS_vec[1] = deltaS
#
#         return hcat(row_dat, row_js, dS_vec)
#     end
#
#     M = grid_map(f, (mx*my,3), grd)
#     dS = M[1:mx*my:end,3]
#     dat = M[:,1]
#     js = Int64[M[:,2]...]
#     is = vcat([k*ones(Int64, mx*my) for k=1:grd.Nk]...)
#
#     matrix = sparse(is, js, dat, grd.Nk, grd._Nx*grd._Ny) * transpose(grd.Proj)
#     return matrix, dS
# end
#
#
# function fieldline_derivatives(bx::Function, by::Function, bz::Function, ds::Float64, mx::Int64, my::Int64, MinvT::Matrix{Float64}, Nz::Int64, grd::Grid)
#
#     # compute forward and backward maps
#     FM, dS1 = fieldline_map_matrix(bx, by, bz, ds, mx, my, MinvT, Nz, grd)
#     @info "Calculated forward map."
#     BM, dS2 = fieldline_map_matrix(bx, by, bz, -ds, mx, my, MinvT, Nz, grd)
#     @info "Calculated backward map."
#
#     # make diagonal matrices for distance coefficients
#     _A = Diagonal(-dS1 ./ (dS2.^2 + dS1.*dS2))
#     _B = Diagonal(1 ./ dS2 - 1 ./ dS1)
#     _C = Diagonal(dS2 ./ (dS1.^2 + dS1.*dS2))
#     _D = Diagonal(1 ./ (dS2.^2 + dS1.*dS2))
#     _E = Diagonal(-1 ./ (dS1.*dS2))
#     _F = Diagonal(1 ./ (dS1.^2 + dS1.*dS2))
#
#     # make off diagonal matrices for kronecker product
#     Ia = spdiagm(-1 => ones(Float64, Nz-1))
#     Ia[1,end] = 1
#     Ib = sparse(I, Nz, Nz)
#     Ic = spdiagm(1 => ones(Float64, Nz-1))
#     Ic[end,1] = 1
#
#     # construct derivative matrices
#     Ds = kron(Ia, _A*BM) + kron(Ib, _B) + kron(Ic, _C*FM)
#     Dss = kron(Ia, _D*BM) + kron(Ib, _E) + kron(Ic, _F*FM)
#
#     return Ds, Dss
# end


function flux_surface_average(psi::Function, psi_val, Ntheta, r0::Vector{Float64}, mx::Int64, my::Int64, MinvT::Matrix{Float64}, grd::Grid; dr=0.01)

    points = zeros(2, Ntheta)
    thetas = LinRange(0, 2*pi, Ntheta+1)
    for i=1:Ntheta

        et = Float64[cos(thetas[i]), sin(thetas[i])]
        r = r0 + dr*et
        sgn = sign(psi(r...) - psi_val)

        while true
            y = r + dr*et
            if sgn*(psi(y...) - psi_val) < 0
                points[:,i] = 0.5*(y + r)
                break
            end
            r = y
        end

    end

    return line_average(points, mx, my, MinvT, grd)
end
