
struct Assets

    GRID :: Grid
    Dx :: SparseMatrixCSC
    Dy :: SparseMatrixCSC
    Dxy :: SparseMatrixCSC
    Dxx :: SparseMatrixCSC
    Dyy :: SparseMatrixCSC
    Dxxx :: SparseMatrixCSC
    Dyyy :: SparseMatrixCSC
    Dxxy :: SparseMatrixCSC
    Dxyy :: SparseMatrixCSC
    Ds :: SparseMatrixCSC
    Dss :: SparseMatrixCSC
    DIFF_lnn :: LinearLeftHandSide
    DIFF_lnTe :: LinearLeftHandSide
    DIFF_lnTi :: LinearLeftHandSide
    DIFF_u :: LinearLeftHandSide
    DIFF_w :: LinearLeftHandSide
    DIFF_A :: LinearLeftHandSide
    HHOLTZ :: LinearLeftHandSide
    P0 :: SparseMatrixCSC
    P1 :: SparseMatrixCSC
    P2 :: SparseMatrixCSC
    P3 :: SparseMatrixCSC
    R1 :: SparseMatrixCSC
    R2 :: SparseMatrixCSC
    R3 :: SparseMatrixCSC
    LAM :: SparseMatrixCSC
    DCHLT1 :: SparseMatrixCSC
    DCHLT2 :: SparseMatrixCSC
    DCHLT3 :: SparseMatrixCSC
    NMANN1 :: SparseMatrixCSC
    NMANN2 :: SparseMatrixCSC
    NMANN3 :: SparseMatrixCSC
    FLXAVG :: SparseMatrixCSC
    TRGT :: Vector{Float64}
    Sn :: Vector{Float64}
    STe :: Vector{Float64}
    STi :: Vector{Float64}
    params :: Vector{Float64}
    dt :: Float64
    N_subcycle :: Int64

end


function Assets(geoinit::Function, phys_ref::Vector{Float64}, Nz::Int64, ds::Float64, dt::Float64; k=0.2, w=0.66, mx=5, my=5, N_subcycle=100)

    # return geoinit()
    grd, bdry, lcfs_pts, trgt_splitter, srcs, bx, by, bz = geoinit()
    params = dimensionless_parameters(phys_ref...)

    P0, P1, P2, P3 = bdry[1]
    R1, R2, R3 = bdry[2]
    DCHLT1, DCHLT2, DCHLT3 = bdry[3]
    NMANN1, NMANN2, NMANN3 = bdry[4]
    Sn, STe, STi = srcs

    # generate stencils
    Mxinv = stencil1d(mx)
    Myinv = stencil1d(my)
    MxyinvT = stencil2d(mx, my)

    # generate derivatives
    Dx = derivative_matrix(1, 0, Mxinv, Myinv, grd)
    Dy = derivative_matrix(0, 1, Mxinv, Myinv, grd)
    Dxy = derivative_matrix(1, 1, Mxinv, Myinv, grd)
    Dxx = derivative_matrix(2, 0, Mxinv, Myinv, grd)
    Dyy = derivative_matrix(0, 2, Mxinv, Myinv, grd)
    Dxxx = derivative_matrix(3, 0, Mxinv, Myinv, grd)
    Dyyy = derivative_matrix(0, 3, Mxinv, Myinv, grd)
    Dxxy = derivative_matrix(2, 1, Mxinv, Myinv, grd)
    Dxyy = derivative_matrix(1, 2, Mxinv, Myinv, grd)
    Ds, Dss = fieldline_derivatives(bx, by, bz, ds, mx, my, MxyinvT, Nz, grd)

    @info "Generated derivative matrices."

    # compute LAM
    LAM = Diagonal(1 ./ (1 .+ dt*diag(P1 + P2 + P3)))

    # setup linear problems
    de2 = params[8]
    Ldiff = I - k*(Dxx + Dyy)
    Lhelm = I - de2*(Dxx + Dyy)
    DIFF_lnn = LinearLeftHandSide(P0*Ldiff + NMANN1 + NMANN2 + NMANN3, w)
    DIFF_lnTe = LinearLeftHandSide(P0*Ldiff + NMANN1 + NMANN2 + NMANN3, w)
    DIFF_lnTi = LinearLeftHandSide(P0*Ldiff + NMANN1 + NMANN2 + NMANN3, w)
    DIFF_u = LinearLeftHandSide(P0*Ldiff + NMANN1 + NMANN2 + DCHLT3, w)
    DIFF_w = LinearLeftHandSide(P0*Ldiff + DCHLT1 + DCHLT2 + DCHLT3, w)
    DIFF_A = LinearLeftHandSide(P0*Ldiff + DCHLT1 + DCHLT2 + DCHLT3, w)
    HHOLTZ = LinearLeftHandSide(P0*Lhelm + DCHLT1 + DCHLT2 + DCHLT3, w)

    @info "Setup linear problems."

    # compute vector for averaging over the LCFS
    rows = SparseMatrixCSC[]
    for i=1:size(lcfs_pts)[2]
        dat, js = interpolation_row(lcfs_pts[:,i]..., mx, my, MxyinvT, grd)
        is = ones(Int64, length(js))
        row = sparse(is, js, dat, 1, grd._Nx*grd._Ny) * transpose(grd.Proj)
        append!(rows, [row])
    end
    FLXAVG = +(rows...) / size(lcfs_pts)[2]

    # define TRGT
    TRGT = f_to_grid((x, y) -> sign(poly_dist(x, y, trgt_splitter)), grd)

    return Assets(grd, Dx, Dy, Dxy, Dxx, Dyy, Dxxx, Dyyy, Dxxy, Dxyy, Ds, Dss,
                     DIFF_lnn, DIFF_lnTe, DIFF_lnTi, DIFF_u, DIFF_w, DIFF_A,
                     HHOLTZ, P0, P1, P2, P3, R1, R2, R3, LAM, DCHLT1, DCHLT2,
                     DCHLT3, NMANN1, NMANN2, NMANN3, FLXAVG, TRGT, Sn, STe,
                     STi, params, dt, N_subcycle)
end


function simple_solovev_assets(Nx, Ny, A, ep, ka, del, L, R0, n0, T0, B0)

    phys_ref = [ep*R0, n0, T0, B0]

    Assets(phys_ref, 1, 0.001, 0.001) do

        psi, bx, by, bz = solovev_flux_function(A, del, ep, ka, 1.0)
        psi_in = psi(1 - ep + L, 0)
        psi_out = 0

        inner_flux_surface = Barrier() do
            rmap(x, y) = trace_reflection(x, y, psi, bx, by, psi_in, 0.001)
            psi, psi_in, rmap
        end

        outer_flux_surface = Barrier() do
            rmap(x, y) = trace_reflection(x, y, psi, bx, by, psi_out, 0.001)
            psi, psi_out, rmap
        end

        bars = [inner_flux_surface, outer_flux_surface, Everywhere()]
        bigdeltas = [1, -1, 1] * L/2
        grd = Grid(Nx, Ny, 1) do
            r0 = [0.01, -(ka*ep + L)]
            r1 = [1 + ep+L, ka*ep + L]
            (x,y) -> check_if_inside(x, y, bigdeltas, bars), r0, r1
        end

        # return grd

        @info "Generated Grid"

        smalldeltas = bigdeltas / 2.0
        qs = ones(3) * 0.01
        bdry = boundary_operators(3, 3, stencil2d(3, 3), smalldeltas, bars, qs, grd)

        @info "Generated boundary operators."

        lcfs_pts = zeros(2, 30)
        ts = LinRange(0, 2*pi, 31)
        ts = ts[1:end-1]
        alpha = asin(del)
        for i=1:30
            lcfs_pts[1,i] = 1 + (ep-L/2)*cos(ts[i] + alpha*sin(ts[i]))
            lcfs_pts[2,i] = (ep-L/2)*ka*sin(ts[i])
        end

        # trgt_splitter doesn't matter because there is no limiter or diverter here
        trgt_splitter = Float64[0 1; 0 1]

        # define constant source for now
        Sn = f_to_grid((x,y) -> 1/(100*dt), grd)
        STe = f_to_grid((x,y) -> 1/(100*dt), grd)
        STi = f_to_grid((x,y) -> 1/(100*dt), grd)

        grd, bdry, lcfs_pts, trgt_splitter, [Sn, STe, STi], bx, by, bz
    end

end
