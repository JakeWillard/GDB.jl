

function rectangle_assets(Nx, Nz, q, a, R0, k, n0, T0, B0, dt)

    # define flux surfaces
    psi(x, y) = -y
    rmap_in(x, y) = [x, -y]
    rmap_out(x, y) = [x, 2 - y]
    inner_flux = Barrier(psi, 0, rmap_in)
    outer_flux = Barrier(psi, -1, rmap_out)

    # define field-aligned sheaths
    vert = zeros(2, 5)
    vert[:,1] = Float64[0, -10]
    vert[:,2] = Float64[0, 10]
    vert[:,3] = Float64[1, 10]
    vert[:,4] = Float64[1, -10]
    vert[:,5] = Float64[0, -10]
    sheaths = PolyBarrier(vert)

    # define grid
    bigdeltas = Float64[-0.3, 0.3, -0.3]
    bars = Barrier[inner_flux, outer_flux, sheaths]
    r0 = Float64[-0.5, -0.5]
    r1 = Float64[1.5, 1.5]
    grd = Grid((x,y) -> check_if_inside(x, y, bigdeltas, bars), r0, r1, Nx, Nx, Nz)

    @info "Computed Grid" Nk=grd.Nk

    # make 5 points stencils
    m = 5
    Minv = stencil1d(m)
    MinvT = Matrix(stencil2d(m, m))

    # create derivative matrices
    Dx = derivative_matrix(1, 0, Minv, Minv, grd)
    Dy = derivative_matrix(0, 1, Minv, Minv, grd)
    Dxy = derivative_matrix(1, 1, Minv, Minv, grd)
    Dxx = derivative_matrix(2, 0, Minv, Minv, grd)
    Dyy = derivative_matrix(0, 2, Minv, Minv, grd)
    Dxxx = derivative_matrix(3, 0, Minv, Minv, grd)
    Dyyy = derivative_matrix(0, 3, Minv, Minv, grd)
    Dxxy = derivative_matrix(2, 1, Minv, Minv, grd)
    Dxyy = derivative_matrix(1, 2, Minv, Minv, grd)

    @info "Generated x and y derivatives"

    # create field-line derivative matrices
    bx(x, y) = 1 / sqrt(1 + (R0*q/a)^2)
    by(x, y) = 0.0
    bz(x, y) = R0*q*bx(x, y) / a
    Ds, Dss = fieldline_derivatives(bx, by, bz, 0.001, m, m, MinvT, Nz, grd)

    @info "Generated field-line derivatives"

    # generate boundary operators
    smallerdeltas = Float64[-0.1, 0.1, -0.1]
    qs = Float64[0.01, 0.01, 0.01]
    PENs, REFs, DCHLTs, NMANNs = boundary_operators(m, m, MinvT, smallerdeltas, bars, qs, grd)
    P0, P1, P2, P3 = PENs
    R1, R2, R3 = REFs
    DCHLT1, DCHLT2, DCHLT3 = DCHLTs
    NMANN1, NMANN2, NMANN3 = NMANNs

    @info "Generated boundary operators"

    # compute LAM matrix
    LAM = Diagonal(1 ./ (1 .+ dt*diag(P1 + P2 + P3)))

    # setup diffusion problems with diffusion constant k
    L = I - k*(Dxx + Dyy)
    DIFF_lnn = LinearLeftHandSide(P0*L + NMANN1 + NMANN2 + NMANN3, 2/3)
    DIFF_lnTe = LinearLeftHandSide(P0*L + NMANN1 + NMANN2 + NMANN3, 2/3)
    DIFF_lnTi = LinearLeftHandSide(P0*L + NMANN1 + NMANN2 + NMANN3, 2/3)
    DIFF_u = LinearLeftHandSide(P0*L + NMANN1 + NMANN2 + DCHLT3, 2/3)
    DIFF_w = LinearLeftHandSide(P0*L + DCHLT1 + DCHLT2 + DCHLT3, 2/3)
    DIFF_A = LinearLeftHandSide(P0*L + DCHLT1 + DCHLT2 + DCHLT3, 2/3)

    # compute parameters
    params = dimensionless_parameters(a, R0, n0, T0, B0)

    # setup helmoltz problem
    de2 = params[8]
    Lhelm = I - de2*(Dxx + Dyy)
    HHOLTZ = LinearLeftHandSide(P0*Lhelm + DCHLT1 + DCHLT2 + DCHLT3, 2/3)

    @info "Setup linear problems"

    # compute FLXAVG, from middle flux surface
    flxpoints = hcat([Float64[x, 0.5] for x in LinRange(0, 1, 15)]...)
    rows = SparseMatrixCSC[]
    for i=1:15
        dat, js = interpolation_row(flxpoints[:,i]..., m, m, MinvT, grd)
        is = ones(Int64, length(js))
        row = sparse(is, js, dat, 1, grd._Nx*grd._Ny) * transpose(grd.Proj)
        append!(rows, [row])
    end
    FLXAVG = +(rows...) / 15 # XXX

    # define TRGT
    TRGT = f_to_grid((x, y) -> sign(x - 0.5), grd)

    # define source terms
    Sn = f_to_grid((x,y) -> 1/(100*dt), grd)
    STe = f_to_grid((x,y) -> 1/(100*dt), grd)
    STi = f_to_grid((x,y) -> 1/(100*dt), grd)

    # subcycle 20 times per timestep
    N_subcycle = 20

    # placeholder seed for now
    seed = Float64[0, 0, 0]

    return Assets(grd, Dx, Dy, Dxy, Dxx, Dyy, Dxxx, Dyyy, Dxxy, Dxyy, Ds, Dss,
                     DIFF_lnn, DIFF_lnTe, DIFF_lnTi, DIFF_u, DIFF_w, DIFF_A,
                     HHOLTZ, P0, P1, P2, P3, R1, R2, R3, LAM, DCHLT1, DCHLT2,
                     DCHLT3, NMANN1, NMANN2, NMANN3, FLXAVG, TRGT, Sn, STe,
                     STi, params, dt, N_subcycle, seed)
end
