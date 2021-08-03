
# store psis in vector and iterate to do derivatives

function solovev_psi_functions()
    # up down symmetric homogeneous solutions
    psi1(x,y) = 1
    psi2(x,y) = x^2
    psi3(x,y) = y^2 - x^2*log(x)
    psi4(x,y) = x^4 - 4*x^2*y^2
    psi5(x,y) = 2*y^4 - 9*y^2*x^2 + 3*x^4*log(x) - 12*x^2*y^2*log(x)
    psi6(x,y) = x^6 - 12*x^4*y^2 + 8*x^2*y^4
    psi7(x,y) = 8*y^6 - 140*y^4*x^2 + 75*y^2*x^4 - 15*x^6*log(x) + 180*x^4*y^2*log(x) - 120*x^2*y^4*log(x)

    # up down asymmetric homogenous solutions
    psi8(x,y) = y
    psi9(x,y) = y*x^2
    psi10(x,y) = y^3 - 3*y*x^2*log(x)
    psi11(x,y) = 3*y*x^4 - 4*y^3*x^2
    psi12(x,y) = 8*y^5 - 45*y*x^4 - 80*y^3*x^2*log(x) + 60*y*x^4*log(x)

    psibasis = [psi1 psi2 psi3 psi4 psi5 psi6 psi7 psi8 psi9 psi10 psi11 psi12]

    return psibasis
end


function solve_solovev_psi_functions(x, y) #XXX might not want to call this function the same as the function above?
    # up down symmetric homogeneous solutions
    psi1 = 1
    psi2 = x^2
    psi3 = y^2 - x^2*log(x)
    psi4 = x^4 - 4*x^2*y^2
    psi5 = 2*y^4 - 9*y^2*x^2 + 3*x^4*log(x) - 12*x^2*y^2*log(x)
    psi6 = x^6 - 12*x^4*y^2 + 8*x^2*y^4
    psi7 = 8*y^6 - 140*y^4*x^2 + 75*y^2*x^4 - 15*x^6*log(x) + 180*x^4*y^2*log(x) - 120*x^2*y^4*log(x)

    # up down asymmetric homogenous solutions
    psi8 = y
    psi9 = y*x^2
    psi10 = y^3 - 3*y*x^2*log(x)
    psi11 = 3*y*x^4 - 4*y^3*x^2
    psi12 = 8*y^5 - 45*y*x^4 - 80*y^3*x^2*log(x) + 60*y*x^4*log(x)

    psivals = [psi1 psi2 psi3 psi4 psi5 psi6 psi7 psi8 psi9 psi10 psi11 psi12]

    return psivals
end


function solovev_parameterize(n, delta, epsilon, kappa)
    # epsilon = reverse aspect ratio (epsilon = a/R0 for R0 the major radius (center of machine to
    # center of region of plasma) and a the minor radius (center of plasma to wall))
    # kappa = elongation
    # sin(alpha) = delta = triangularity
    # n = resolution of points in x/y

    alpha = asin(delta)

    # calculating x and y coordinates of flux surface cross section
    tau = LinRange(0, 2*pi, n)
    x = zeros(n)
    y = zeros(n)
    for i=1:n
        x[i] = 1 + epsilon*cos(tau[i] + alpha*sin(tau[i]))
        y[i] = epsilon*kappa*sin(tau[i])
    end

    return tau, x, y
end


function solovev_particular_psi(A, x, y)
    psi = x^4/8 + A * (1/2*x^2*log(x) - x^4/8)
    return psi
end


function solovev_homo_partial_x(functs, xval, yval) #XXX could be confusing that the input variable psis is also called the same as two seperate functions defined above.
    partial_x = zeros(12)
    for i=1:12
        # partial_x[i] = ForwardDiff.derivative(x -> psis[i](x,y), x)

        # XXX: ForwardDiff.derivative returns a number, not a function. The input looks like (::Function, ::Real), and the function has to be single variable
        #      so it doesn't know what to do with psis[i](x,y), it has to be psis[i](x,yval)
        partial_x[i] = ForwardDiff.derivative(x -> functs[i](x,yval), xval) #XXX <-- applying this correction to the other instances should fix the issue.
    end
    # for j=1:12
    #     partial_x[j] = partial_x[j][xval, yval]
    # end

    return partial_x
end


function solovev_homo_partial_xx(functs, xval, yval)
    # px = []
    pxx = zeros(12)
    for i=1:12
        # XXX Since ForwardDiff isn't returning functions, this is the way to take a second derivative:
        pxx[i] = ForwardDiff.derivative(x -> ForwardDiff.derivative(x -> functs[i](x,yval), x), xval)
    end
    return pxx
end


function solovev_homo_partial_y(functs, xval, yval) #XXX Need to just apply the changes above to the y versions.
    partial_y = zeros(12)
    for i=1:12
        partial_y[i] = ForwardDiff.derivative(y -> functs[i](xval,y), yval)
    end
    return partial_y
end


function solovev_homo_partial_yy(functs, xval, yval)
    # py = zeros(12)
    pyy = zeros(12)
    for i=1:12
        # py[i] = ForwardDiff.derivative(y -> psis[i](x,y), y)
        pyy[i] = ForwardDiff.derivative(y -> ForwardDiff.derivative(y -> functs[i](xval,y), y), yval)
    end
    return pyy
end


function solovev_particular_partial_x(A, xval, yval)
    partial_x = ForwardDiff.derivative(x -> solovev_particular_psi(A, x, yval), xval)
    return partial_x
end


function solovev_particular_partial_xx(A, xval, yval)
    # partial_x = ForwardDiff.derivative(x -> solovev_particular_psi(A, x, y), x)
    partial_xx = ForwardDiff.derivative(x -> ForwardDiff.derivative(x -> solovev_particular_psi(A, x, yval), x), xval)
    return partial_xx
end


function solovev_particular_partial_y(A, xval, yval)
    partial_y = ForwardDiff.derivative(y -> solovev_particular_psi(A, xval, y), yval)
    return partial_y
end


function solovev_particular_partial_yy(A, xval, yval)
    # partial_y = ForwardDiff.derivative(y -> solovev_particular_psi(A)(y), yval)
    partial_yy = ForwardDiff.derivative(y -> ForwardDiff.derivative(y -> solovev_particular_psi(A, xval, y), y), yval)
    return partial_yy
end


function solovev_flux_function(A, delta, epsilon, kappa, B0; upsep=[0,0], downsep=[0,0])
    alpha = asin(delta)

    # last terms are xsep/ysep - can be changed to be inputs if need be
    N1 = -(1 + alpha)^2/(epsilon*kappa^2)
    N2 = (1 - alpha)^2/(epsilon*kappa^2)
    N3 = -kappa/(epsilon*(cos(alpha)^2))

    x = [1 + epsilon, 1 - epsilon, 1 - delta*epsilon, upsep[1], downsep[1]]
    y = [0, 0, kappa*epsilon, upsep[2], downsep[2]]

    psisbasis = solovev_psi_functions()

    rhs1 = solovev_particular_psi(A, x[1], y[1])
    rhsy2 = solovev_particular_partial_y(A, x[1], y[1])
    rhsyy3 = solovev_particular_partial_yy(A, x[1], y[1])
    rhsx3 = solovev_particular_partial_x(A, x[1], y[1])

    rhs4 = solovev_particular_psi(A, x[2], y[2])
    rhsy5 = solovev_particular_partial_y(A, x[2], y[2])
    rhsyy6 = solovev_particular_partial_yy(A, x[2], y[2])
    rhsx6 = solovev_particular_partial_x(A, x[2], y[2])

    psis1 = solve_solovev_psi_functions(x[1], y[1])
    psisy2 = solovev_homo_partial_y(psisbasis, x[1], y[1])
    psisyy3 = solovev_homo_partial_yy(psisbasis, x[1], y[1])
    psisx3 = solovev_homo_partial_x(psisbasis, x[1], y[1])

    psis4 = solve_solovev_psi_functions(x[2], y[2])
    psisy5 = solovev_homo_partial_y(psisbasis, x[2], y[2])
    psisyy6 = solovev_homo_partial_yy(psisbasis, x[2], y[2])
    psisx6 = solovev_homo_partial_x(psisbasis, x[2], y[2])

    Amat = zeros(12, 12)
    Amat[1,:] = psis1
    Amat[2,:] = psisy2
    Amat[3,:] = psisyy3 + N1*psisx3
    Amat[4,:] = psis4
    Amat[5,:] = psisy5
    Amat[6,:] = psisyy6 + N2*psisx6

    if upsep==[0,0] && downsep==[0,0]
        # calculating values of flux function and derivatives at boundaries
        rhs7 = solovev_particular_psi(A, x[3], y[3])
        rhsx8 = solovev_particular_partial_x(A, x[3], y[3])
        rhsxx9 = solovev_particular_partial_xx(A, x[3], y[3])
        rhsy9 = solovev_particular_partial_y(A, x[3], y[3])
        rhs10 = solovev_particular_psi(A, x[3], -y[3])
        rhsx11 = solovev_particular_partial_x(A, x[3], -y[3])
        rhsxx12 = solovev_particular_partial_xx(A, x[3], -y[3])
        rhsy12 = solovev_particular_partial_y(A, x[3], -y[3])

        psis7 = solve_solovev_psi_functions(x[3], y[3])
        psisx8 = solovev_homo_partial_x(psisbasis, x[3], y[3])
        psisxx9 = solovev_homo_partial_xx(psisbasis, x[3], y[3])
        psisy9 = solovev_homo_partial_y(psisbasis, x[3], y[3])
        psis10 = solve_solovev_psi_functions(x[3], -y[3])
        psisx11 = solovev_homo_partial_x(psisbasis, x[3], -y[3])
        psisxx12 = solovev_homo_partial_xx(psisbasis, x[3], -y[3])
        psisy12 = solovev_homo_partial_y(psisbasis, x[3], -y[3])

        Amat[7,:] = psis7
        Amat[8,:] = psisx8
        Amat[9,:] = psisxx9 + N3*psisy9
        Amat[10,:] = psis10
        Amat[11,:] = psisx11
        Amat[12,:] = psisxx12 - N3*psisy12

        b = [rhs1, rhsy2, rhsyy3+N1*rhsx3, rhs4, rhsy5, rhsyy6+N2*rhsx6,
        rhs7, rhsx8, rhsxx9+N3*rhsy9, rhs10, rhsx11, rhsxx12-N3*rhsy12]
    elseif downsep==[0,0]
        rhs7 = solovev_particular_psi(A, x[4], y[4])
        rhsx8 = solovev_particular_partial_x(A, x[4], y[4])
        rhsy9 = solovev_particular_partial_y(A, x[4], y[4])
        rhs10 = solovev_particular_psi(A, x[3], -y[3])
        rhsx11 = solovev_particular_partial_x(A, x[3], -y[3])
        rhsxx12 = solovev_particular_partial_xx(A, x[3], -y[3])
        rhsy12 = solovev_particular_partial_y(A, x[3], -y[3])

        psis7 = solve_solovev_psi_functions(x[4], y[4])
        psisx8 = solovev_homo_partial_x(psisbasis, x[4], y[4])
        psisy9 = solovev_homo_partial_y(psisbasis, x[4], y[4])
        psis10 = solve_solovev_psi_functions(x[3], -y[3])
        psisx11 = solovev_homo_partial_x(psisbasis, x[3], -y[3])
        psisxx12 = solovev_homo_partial_xx(psisbasis, x[3], -y[3])
        psisy12 = solovev_homo_partial_y(psisbasis, x[3], -y[3])

        Amat[7,:] = psis7
        Amat[8,:] = psisx8
        Amat[9,:] = psisy9
        Amat[10,:] = psis10
        Amat[11,:] = psisx11
        Amat[12,:] = psisxx12 - N3*psisy12

        b = [rhs1, rhsy2, rhsyy3+N1*rhsx3, rhs4, rhsy5, rhsyy6+N2*rhsx6,
        rhs7, rhsx8, rhsy9, rhs10, rhsx11, rhsxx12+N3*rhsy12]
    elseif upsep==[0,0]
        rhs7 = solovev_particular_psi(A, x[3], y[3])
        rhsx8 = solovev_particular_partial_x(A, x[3], y[3])
        rhsxx9 = solovev_particular_partial_xx(A, x[3], y[3])
        rhsy9 = solovev_particular_partial_y(A, x[3], y[3])
        rhs10 = solovev_particular_psi(A, x[5], y[5])
        rhsx11 = solovev_particular_partial_x(A, x[5], y[5])
        rhsy12 = solovev_particular_partial_y(A, x[5], y[5])

        psis7 = solve_solovev_psi_functions(x[3], y[3])
        psisx8 = solovev_homo_partial_x(psisbasis, x[3], y[3])
        psisxx9 = solovev_homo_partial_xx(psisbasis, x[3], y[3])
        psisy9 = solovev_homo_partial_y(psisbasis, x[3], y[3])
        psis10 = solve_solovev_psi_functions(x[5], y[5])
        psisx11 = solovev_homo_partial_x(psisbasis, x[5], y[5])
        psisy12 = solovev_homo_partial_y(psisbasis, x[5], y[5])

        Amat[7,:] = psis7
        Amat[8,:] = psisx8
        Amat[9,:] = psisxx9 + N3*psisy9
        Amat[10,:] = psis10
        Amat[11,:] = psisx11
        Amat[12,:] = psisy12

        b = [rhs1, rhsy2, rhsyy3+N1*rhsx3, rhs4, rhsy5, rhsyy6+N2*rhsx6,
        rhs7, rhsx8, rhsxx9+N3*rhsy9, rhs10, rhsx11, rhsy12]
    else
        rhs7 = solovev_particular_psi(A, x[4], y[4])
        rhsx8 = solovev_particular_partial_x(A, x[4], y[4])
        rhsy9 = solovev_particular_partial_y(A, x[4], y[4])
        rhs10 = solovev_particular_psi(A, x[5], y[5])
        rhsx11 = solovev_particular_partial_x(A, x[5], y[5])
        rhsy12 = solovev_particular_partial_y(A, x[5], y[5])

        psis7 = solve_solovev_psi_functions(x[4], y[4])
        psisx8 = solovev_homo_partial_x(psisbasis, x[4], y[4])
        psisy9 = solovev_homo_partial_y(psisbasis, x[4], y[4])
        psis10 = solve_solovev_psi_functions(x[5], y[5])
        psisx11 = solovev_homo_partial_x(psisbasis, x[5], y[5])
        psisy12 = solovev_homo_partial_y(psisbasis, x[5], y[5])

        Amat[7,:] = psis7
        Amat[8,:] = psisx8
        Amat[9,:] = psisy9
        Amat[10,:] = psis10
        Amat[11,:] = psisx11
        Amat[12,:] = psisy12

        b = [rhs1, rhsy2, rhsyy3+N1*rhsx3, rhs4, rhsy5, rhsyy6+N2*rhsx6,
        rhs7, rhsx8, rhsy9, rhs10, rhsx11, rhsy12]
    end

    c = Amat \ -b
    psifunctions = solovev_psi_functions()
    psi(x,y) = solovev_particular_psi(A,x,y) .+ dot(c, [p(x,y) for p in psifunctions])

    Bx(x, y) = ForwardDiff.derivative(u -> -psi(x,u), y)
    By(x, y) = ForwardDiff.derivative(u -> psi(u,y), x) #XXX changed the sign of this to make consistent with trace.jl, should check later which way is actually correct. 
    # B0 = R0^4*B0^2/psi0^2
    Bphi(x,y) = (B0^2 - 2*A*psi(x,y))^(1/2)
    B(x,y) = norm(Float64[Bx(x,y), By(x,y), Bphi(x,y)])

    bx(x,y) = Bx(x,y) / B(x,y)
    by(x,y) = By(x,y) / B(x,y)
    bz(x,y) = Bpsi(x,y) / B(x,y)

    return psi, bx, by, bz
end
