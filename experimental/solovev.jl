using Base: Float64
# solutions to Grad-Shafranov using Solov'ev profiles
using DifferentialEquations
using Plots
using ForwardDiff
using LinearAlgebra

# store phis in vector and iterate to do derivatives

function phifuncts()
    # up down symmetric homogeneous solutions


    phi1(x,y) = 1
    phi2(x,y) = x^2
    phi3(x,y) = y^2 - x^2*log(x)
    phi4(x,y) = x^4 - 4*x^2*y^2
    phi5(x,y) = 2*y^4 - 9*y^2*x^2 + 3*x^4*log(x) - 12*x^2*y^2*log(x)
    phi6(x,y) = x^6 - 12*x^4*y^2 + 8*x^2*y^4
    phi7(x,y) = 8*y^6 - 140*y^4*x^2 + 75*y^2*x^4 - 15*x^6*log(x) + 180*x^4*y^2*log(x) - 120*x^2*y^4*log(x)

    # up down asymmetric homogenous solutions
    phi8(x,y) = y
    phi9(x,y) = y*x^2
    phi10(x,y) = y^3 - 3*y*x^2*log(x)
    phi11(x,y) = 3*y*x^4 - 4*y^3*x^2
    phi12(x,y) = 8*y^5 - 45*y*x^4 - 80*y^3*x^2*log(x) + 60*y*x^4*log(x)

    phibasis = [phi1 phi2 phi3 phi4 phi5 phi6 phi7 phi8 phi9 phi10 phi11 phi12]

    return phibasis
end


function phis(x, y) #XXX might not want to call this function the same as the function above?
    # up down symmetric homogeneous solutions
    phi1 = 1
    phi2 = x^2
    phi3 = y^2 - x^2*log(x)
    phi4 = x^4 - 4*x^2*y^2
    phi5 = 2*y^4 - 9*y^2*x^2 + 3*x^4*log(x) - 12*x^2*y^2*log(x)
    phi6 = x^6 - 12*x^4*y^2 + 8*x^2*y^4
    phi7 = 8*y^6 - 140*y^4*x^2 + 75*y^2*x^4 - 15*x^6*log(x) + 180*x^4*y^2*log(x) - 120*x^2*y^4*log(x)

    # up down asymmetric homogenous solutions
    phi8 = y
    phi9 = y*x^2
    phi10 = y^3 - 3*y*x^2*log(x)
    phi11 = 3*y*x^4 - 4*y^3*x^2
    phi12 = 8*y^5 - 45*y*x^4 - 80*y^3*x^2*log(x) + 60*y*x^4*log(x)

    phivals = [phi1 phi2 phi3 phi4 phi5 phi6 phi7 phi8 phi9 phi10 phi11 phi12]

    return phivals
end


function parameterize(n, delta, epsilon, kappa)
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

function parameterize_extra(n, m, delta, epsilon, kappa)
    alpha = asin(delta)

    # calculating x and y coordinates of flux surface cross section
    tau = LinRange(0, 2*pi, n)
    epsres = LinRange(0, epsilon, m)
    x = zeros(n*m)
    y = zeros(n*m)
    for i=1:n
        for j=1:m
            x[(i-1)*m+j] = 1 + epsres[j]*cos(tau[i] + alpha*sin(tau[i]))
            y[(i-1)*m+j] = epsres[j]*kappa*sin(tau[i])
            # x[i] = 1 + epsilon*cos(tau[i] + alpha*sin(tau[i]))
            # y[i] = epsilon*kappa*sin(tau[i])
        end

    end

    return tau, x, y
end


function calcphi(A, x, y)
    phi = x^4/8 + A * (1/2*x^2*log(x) - x^4/8)
    return phi
end


function calcphifunct(A)
    phi(x,y) = x^4/8 + A * (1/2*x^2*log(x) - x^4/8)
    return phi
end


function phisdiffx(functs, xval, yval) #XXX could be confusing that the input variable phis is also called the same as two seperate functions defined above.
    partial_x = zeros(12)
    for i=1:12
        # partial_x[i] = ForwardDiff.derivative(x -> phis[i](x,y), x)

        # XXX: ForwardDiff.derivative returns a number, not a function. The input looks like (::Function, ::Real), and the function has to be single variable
        #      so it doesn't know what to do with phis[i](x,y), it has to be phis[i](x,yval)
        partial_x[i] = ForwardDiff.derivative(x -> functs[i](x,yval), xval) #XXX <-- applying this correction to the other instances should fix the issue.
    end
    # for j=1:12
    #     partial_x[j] = partial_x[j][xval, yval]
    # end

    return partial_x
end


function phisdiffxx(functs, xval, yval)
    # px = []
    pxx = zeros(12)
    for i=1:12
        # XXX Since ForwardDiff isn't returning functions, this is the way to take a second derivative:
        pxx[i] = ForwardDiff.derivative(x -> ForwardDiff.derivative(x -> functs[i](x,yval), x), xval)
    end
    return pxx
end


function phisdiffy(functs, xval, yval) #XXX Need to just apply the changes above to the y versions.
    partial_y = zeros(12)
    for i=1:12
        partial_y[i] = ForwardDiff.derivative(y -> functs[i](xval,y), yval)
    end
    return partial_y
end


function phisdiffyy(functs, xval, yval)
    # py = zeros(12)
    pyy = zeros(12)
    for i=1:12
        # py[i] = ForwardDiff.derivative(y -> phis[i](x,y), y)
        pyy[i] = ForwardDiff.derivative(y -> ForwardDiff.derivative(y -> functs[i](xval,y), y), yval)
    end
    return pyy
end


function diffx(A, xval, yval)
    partial_x = ForwardDiff.derivative(x -> calcphifunct(A)(x, yval), xval)
    return partial_x
end


function diffxx(A, xval, yval)
    # partial_x = ForwardDiff.derivative(x -> calcphi(A, x, y), x)
    partial_xx = ForwardDiff.derivative(x -> ForwardDiff.derivative(x -> calcphifunct(A)(x, yval), x), xval)
    return partial_xx
end


function diffy(A, xval, yval)
    partial_y = ForwardDiff.derivative(y -> calcphifunct(A)(xval, y), yval)
    return partial_y
end


function diffyy(A, xval, yval)
    # partial_y = ForwardDiff.derivative(y -> calcphi(A)(y), yval)
    partial_yy = ForwardDiff.derivative(y -> ForwardDiff.derivative(y -> calcphifunct(A)(xval, y), y), yval)
    return partial_yy
end


function solve_boundaries(A, delta, epsilon, kappa)
    # complete sets of boundary values at which we solve for values of c's
    alpha = asin(delta)
    x = [1 + epsilon, 1 - epsilon, 1 - delta*epsilon]
    y = [0, 0, kappa*epsilon]

    N1 = -(1 + alpha)^2/(epsilon*kappa^2)
    N2 = (1 - alpha)^2/(epsilon*kappa^2)
    N3 = -kappa/(epsilon*(cos(alpha)^2))

    # calculating values of flux function and derivatives at boundaries
    rhsphi1 = calcphi(A, x[1], y[1])
    rhsphi2 = calcphi(A, x[2], y[2])
    rhsphi3 = calcphi(A, x[3], y[3])
    rhsx1 = diffx(A, x[3], y[3])
    rhsyy2 = diffyy(A, x[1], y[1])
    rhsx2 = diffx(A, x[1], y[1])
    rhsyy3 = diffyy(A, x[2], y[2])
    rhsx3 = diffx(A, x[2], y[2])
    rhsxx4 = diffxx(A, x[3], y[3])
    rhsy4 = diffy(A, x[3], y[3])

    phisbasis = phifuncts()

    phis1 = phis(x[1], y[1])
    phis2 = phis(x[2], y[2])
    phis3 = phis(x[3], y[3])
    phisx1 = phisdiffx(phisbasis, x[3], y[3])
    phisyy2 = phisdiffyy(phisbasis, x[1], y[1])
    phisx2 = phisdiffx(phisbasis, x[1], y[1])
    phisyy3 = phisdiffyy(phisbasis, x[2], y[2])
    phisx3 = phisdiffx(phisbasis, x[2], y[2])
    phisxx4 = phisdiffxx(phisbasis, x[3], y[3])
    phisy4 = phisdiffy(phisbasis, x[3], y[3])


    # constructing matrices to solve equation A*c = b
    A = zeros(7,7)
    A[1,:] = phis1[1:7]
    A[2,:] = phis2[1:7]
    A[3,:] = phis3[1:7]
    A[4,:] = phisx1[1:7]
    A[5,:] = phisyy2[1:7] + N1*phisx2[1:7]
    A[6,:] = phisyy3[1:7] + N2*phisx3[1:7]
    A[7,:] = phisxx4[1:7] + N3*phisy4[1:7]

    # need to multiply b by -1 since it's changing sides
    b = [rhsphi1, rhsphi2, rhsphi3, rhsx1, rhsyy2+N1*rhsx2, rhsyy3+N2*rhsx3, rhsxx4+N3*rhsy4]

    c = A \ -b

    return c
end


function boundaries_separatrix(A, delta, epsilon, kappa)
    # doing the same as above with the separatrix boundary conditions
    alpha = asin(delta)
    x = [1 + epsilon, 1 - epsilon, 1 - 1.1*delta*epsilon]
    y = [0, 0, 1.1*kappa*epsilon]

    N1 = -(1 + alpha)^2/(epsilon*kappa^2)
    N2 = (1 - alpha)^2/(epsilon*kappa^2)

    rhsphi1 = calcphi(A, x[1], y[1])
    rhsphi2 = calcphi(A, x[2], y[2])
    rhsphi3 = calcphi(A, x[3], y[3])
    rhsx1 = diffx(A, x[3], y[3])
    rhsy2 = diffy(A, x[3], y[3])
    rhsyy3 = diffyy(A, x[1], y[1])
    rhsx3 = diffx(A, x[1], y[1])
    rhsyy4 = diffyy(A, x[2], y[2])
    rhsx4 = diffx(A, x[2], y[2])

    phisbasis = phifuncts()

    phis1 = phis(x[1], y[1])
    phis2 = phis(x[2], y[2])
    phis3 = phis(x[3], y[3])
    phisx1 = phisdiffx(phisbasis, x[3], y[3])
    phisy2 = phisdiffy(phisbasis, x[3], y[3])
    phisyy3 = phisdiffyy(phisbasis, x[1], y[1])
    phisx3 = phisdiffx(phisbasis, x[1], y[1])
    phisyy4 = phisdiffyy(phisbasis, x[2], y[2])
    phisx4 = phisdiffx(phisbasis, x[2], y[2])

    A = zeros(7,7)
    A[1,:] = phis1[1:7]
    A[2,:] = phis2[1:7]
    A[3,:] = phis3[1:7]
    A[4,:] = phisx1[1:7]
    A[5,:] = phisy2[1:7]
    A[6,:] = phisyy3[1:7] + N1*phisx3[1:7]
    A[7,:] = phisyy4[1:7] + N2*phisx4[1:7]

    b = [rhsphi1, rhsphi2, rhsphi3, rhsx1, rhsy2, rhsyy3+N1*rhsx3, rhsyy4+N2*rhsx4]

    c = A \ -b
    return c
end


function asymmetric_boundaries(A, delta, epsilon, kappa)
    alpha = asin(delta)

    # last terms are xsep/ysep - can be changed to be inputs if need be
    x = [1 + epsilon, 1 - epsilon, 1 - delta*epsilon, 1 - 1.1*delta*epsilon]
    y = [0, 0, kappa*epsilon, -1.1*kappa*epsilon]

    N1 = -(1 + alpha)^2/(epsilon*kappa^2)
    N2 = (1 - alpha)^2/(epsilon*kappa^2)
    N3 = -kappa/(epsilon*(cos(alpha)^2))

    rhsphi1 = calcphi(A, x[1], y[1])
    rhsphi2 = calcphi(A, x[2], y[2])
    rhsphi3 = calcphi(A, x[3], y[3])
    rhsphi4 = calcphi(A, x[4], y[4])
    rhsy5 = diffy(A, x[1], y[1])
    rhsy6 = diffy(A, x[2], y[2])
    rhsx7 = diffx(A, x[3], y[3])
    rhsx8 = diffx(A, x[4], y[4]) # B_y = 0 at lower X-point
    rhsy9 = diffy(A, x[4], y[4]) # B_x = 0 at lower X-point
    rhsyy10 = diffyy(A, x[1], y[1])
    rhsx10 = diffx(A, x[1], y[1])
    rhsyy11 = diffyy(A, x[2], y[2])
    rhsx11 = diffx(A, x[2], y[2])
    rhsxx12 = diffxx(A, x[3], y[3])
    rhsy12 = diffy(A, x[3], y[3])

    phisbasis = phifuncts()

    phis1 = phis(x[1], y[1])
    phis2 = phis(x[2], y[2])
    phis3 = phis(x[3], y[3])
    phis4 = phis(x[4], y[4])
    phisy5 = phisdiffy(phisbasis, x[1], y[1])
    phisy6 = phisdiffy(phisbasis, x[2], y[2])
    phisx7 = phisdiffx(phisbasis, x[3], y[3])
    phisx8 = phisdiffx(phisbasis, x[4], y[4])
    phisy9 = phisdiffy(phisbasis, x[4], y[4])
    phisyy10 = phisdiffyy(phisbasis, x[1], y[1])
    phisx10 = phisdiffx(phisbasis, x[1], y[1])
    phisyy11 = phisdiffyy(phisbasis, x[2], y[2])
    phisx11 = phisdiffx(phisbasis, x[2], y[2])
    phisxx12 = phisdiffxx(phisbasis, x[3], y[3])
    phisy12 = phisdiffy(phisbasis, x[3], y[3])

    A = zeros(12, 12)
    A[1,:] = phis1
    A[2,:] = phis2
    A[3,:] = phis3
    A[4,:] = phis4
    A[5,:] = phisy5
    A[6,:] = phisy6
    A[7,:] = phisx7
    A[8,:] = phisx8
    A[9,:] = phisy9
    A[10,:] = phisyy10 + N1*phisx10
    A[11,:] = phisyy11 + N2*phisx11
    A[12,:] = phisxx12 + N3*phisy12

    b = [rhsphi1, rhsphi2, rhsphi3, rhsphi4, rhsy5, rhsy6, rhsx7, rhsx8, rhsy9,
    rhsyy10+N1*rhsx10, rhsyy11+N2*rhsx11, rhsxx12+N3*rhsy12]

    c = A \ -b

    return c
end


function circle(A)
    radians = LinRange(0, 2*pi*10/11, 11)
    x = zeros(size(radians))
    y = zeros(size(radians))
    for i=1:11
        x[i] = 1 + .5*cos(radians[i])
        y[i] = .5*sin(radians[i])
    end

    rhs1 = calcphi(A, x[1], y[1])
    rhs2 = calcphi(A, x[2], y[2])
    rhs3 = calcphi(A, x[3], y[3])
    rhs4 = calcphi(A, x[4], y[4])
    rhs5 = calcphi(A, x[5], y[5])
    rhs6 = calcphi(A, x[6], y[6])
    rhs7 = calcphi(A, x[7], y[7])
    rhs8 = calcphi(A, x[8], y[8])
    rhs9 = calcphi(A, x[9], y[9])
    rhs10 = calcphi(A, x[10], y[10])
    rhs11 = calcphi(A, x[11], y[11])
    rhsy1 = diffy(A, x[1], y[1])

    phisbasis = phifuncts()

    phis1 = phis(x[1], y[1])
    phis2 = phis(x[2], y[2])
    phis3 = phis(x[3], y[3])
    phis4 = phis(x[4], y[4])
    phis5 = phis(x[5], y[5])
    phis6 = phis(x[6], y[6])
    phis7 = phis(x[7], y[7])
    phis8 = phis(x[8], y[8])
    phis9 = phis(x[9], y[9])
    phis10 = phis(x[10], y[10])
    phis11 = phis(x[11], y[11])
    phisy1 = phisdiffy(phisbasis, x[1], y[1])

    A = zeros(12, 12)
    A[1,:] = phis1
    A[2,:] = phis2
    A[3,:] = phis3
    A[4,:] = phis4
    A[5,:] = phis5
    A[6,:] = phis6
    A[7,:] = phis7
    A[8,:] = phis8
    A[9,:] = phis9
    A[10,:] = phis10
    A[11,:] = phis11
    A[12,:] = phisy1

    b = [rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7, rhs8, rhs9, rhs10, rhs11, rhsy1]

    c = A \ -b
    plot()

    for dummy=1:11
        x[dummy] = x[dummy]*100
        y[dummy] = (y[dummy] + 1)*100
    end
    scatter!(x,y)


    return c
end


function generate_equation(A, delta, epsilon, kappa)
    cvals = solve_boundaries(A, delta, epsilon, kappa)
    # println(cvals)
    phi = phifuncts()
    particularphi = calcphifunct(A)
    # print(phi)

    f(x,y) = particularphi(x,y) .+ dot(cvals[1:7], [p(x,y) for p in phi[1:7]])

    return f
end


function gen_separatrix(A, delta, epsilon, kappa)
    cvals = boundaries_separatrix(A, delta, epsilon, kappa)
    phi = phifuncts()
    particularphi = calcphifunct(A)

    f(x,y) = particularphi(x,y) .+ dot(cvals[1:7], [p(x,y) for p in phi[1:7]])

    return f
end


function gen_asym(A, delta, epsilon, kappa)
    cvals = asymmetric_boundaries(A, delta, epsilon, kappa)
    phi = phifuncts()
    particularphi = calcphifunct(A)

    f(x,y) = particularphi(x,y) .+ dot(cvals, [p(x,y) for p in phi])

    return f
end

function gen_circle(A)
    cvals = circle(A)
    phi = phifuncts()
    particularphi = calcphifunct(A)

    f(x,y) = particularphi(x,y) .+ dot(cvals, [p(x,y) for p in phi])

    return f
end


A = -.155
delta = .33
kappa = 1.7
epsilon = .32


x = LinRange(.25, 1.75, 200)
y = LinRange(-.7, .7, 200)
# tau, x, y = parameterize(200, delta, epsilon, kappa)
z = zeros(Float64, (200, 200))
# psi = gen_asym(A, delta, epsilon, kappa)

# contour(z, levels=1000)
# contour(z, levels = 100, fill=false)
# scatter(x, y)
# plot(x,y, xlim=(0.5, 1.5), ylim = (-.6,.6))
# heatmap(x, y, z)


function gen_flux_function(Aval, delta, epsilon, kappa, upsep=[0,0], downsep=[0,0])
    alpha = asin(delta)

    # last terms are xsep/ysep - can be changed to be inputs if need be
    N1 = -(1 + alpha)^2/(epsilon*kappa^2)
    N2 = (1 - alpha)^2/(epsilon*kappa^2)
    N3 = -kappa/(epsilon*(cos(alpha)^2))

    x = [1 + epsilon, 1 - epsilon, 1 - delta*epsilon, upsep[1], downsep[1]]
    y = [0, 0, kappa*epsilon, upsep[2], downsep[2]]

    phisbasis = phifuncts()

    rhs1 = calcphi(A, x[1], y[1])
    rhsy2 = diffy(A, x[1], y[1])
    rhsyy3 = diffyy(A, x[1], y[1])
    rhsx3 = diffx(A, x[1], y[1])

    rhs4 = calcphi(A, x[2], y[2])
    rhsy5 = diffy(A, x[2], y[2])
    rhsyy6 = diffyy(A, x[2], y[2])
    rhsx6 = diffx(A, x[2], y[2])

    phis1 = phis(x[1], y[1])
    phisy2 = phisdiffy(phisbasis, x[1], y[1])
    phisyy3 = phisdiffyy(phisbasis, x[1], y[1])
    phisx3 = phisdiffx(phisbasis, x[1], y[1])

    phis4 = phis(x[2], y[2])
    phisy5 = phisdiffy(phisbasis, x[2], y[2])
    phisyy6 = phisdiffyy(phisbasis, x[2], y[2])
    phisx6 = phisdiffx(phisbasis, x[2], y[2])

    Amat = zeros(12, 12)
    Amat[1,:] = phis1
    Amat[2,:] = phisy2
    Amat[3,:] = phisyy3 + N1*phisx3
    Amat[4,:] = phis4
    Amat[5,:] = phisy5
    Amat[6,:] = phisyy6 + N2*phisx6

    if upsep==[0,0] && downsep==[0,0]
        # calculating values of flux function and derivatives at boundaries
        rhs7 = calcphi(A, x[3], y[3])
        rhsx8 = diffx(A, x[3], y[3])
        rhsxx9 = diffxx(A, x[3], y[3])
        rhsy9 = diffy(A, x[3], y[3])
        rhs10 = calcphi(A, x[3], -y[3])
        rhsx11 = diffx(A, x[3], -y[3])
        rhsxx12 = diffxx(A, x[3], -y[3])
        rhsy12 = diffy(A, x[3], -y[3])

        phis7 = phis(x[3], y[3])
        phisx8 = phisdiffx(phisbasis, x[3], y[3])
        phisxx9 = phisdiffxx(phisbasis, x[3], y[3])
        phisy9 = phisdiffy(phisbasis, x[3], y[3])
        phis10 = phis(x[3], -y[3])
        phisx11 = phisdiffx(phisbasis, x[3], -y[3])
        phisxx12 = phisdiffxx(phisbasis, x[3], -y[3])
        phisy12 = phisdiffy(phisbasis, x[3], -y[3])

        Amat[7,:] = phis7
        Amat[8,:] = phisx8
        Amat[9,:] = phisxx9 + N3*phisy9
        Amat[10,:] = phis10
        Amat[11,:] = phisx11
        Amat[12,:] = phisxx12 - N3*phisy12

        b = [rhs1, rhsy2, rhsyy3+N1*rhsx3, rhs4, rhsy5, rhsyy6+N2*rhsx6,
        rhs7, rhsx8, rhsxx9+N3*rhsy9, rhs10, rhsx11, rhsxx12-N3*rhsy12]
    elseif downsep==[0,0]
        rhs7 = calcphi(A, x[4], y[4])
        rhsx8 = diffx(A, x[4], y[4])
        rhsy9 = diffy(A, x[4], y[4])
        rhs10 = calcphi(A, x[3], -y[3])
        rhsx11 = diffx(A, x[3], -y[3])
        rhsxx12 = diffxx(A, x[3], -y[3])
        rhsy12 = diffy(A, x[3], -y[3])

        phis7 = phis(x[4], y[4])
        phisx8 = phisdiffx(phisbasis, x[4], y[4])
        phisy9 = phisdiffy(phisbasis, x[4], y[4])
        phis10 = phis(x[3], -y[3])
        phisx11 = phisdiffx(phisbasis, x[3], -y[3])
        phisxx12 = phisdiffxx(phisbasis, x[3], -y[3])
        phisy12 = phisdiffy(phisbasis, x[3], -y[3])

        Amat[7,:] = phis7
        Amat[8,:] = phisx8
        Amat[9,:] = phisy9
        Amat[10,:] = phis10
        Amat[11,:] = phisx11
        Amat[12,:] = phisxx12 - N3*phisy12

        b = [rhs1, rhsy2, rhsyy3+N1*rhsx3, rhs4, rhsy5, rhsyy6+N2*rhsx6,
        rhs7, rhsx8, rhsy9, rhs10, rhsx11, rhsxx12+N3*rhsy12]
    elseif upsep==[0,0]
        rhs7 = calcphi(A, x[3], y[3])
        rhsx8 = diffx(A, x[3], y[3])
        rhsxx9 = diffxx(A, x[3], y[3])
        rhsy9 = diffy(A, x[3], y[3])
        rhs10 = calcphi(A, x[5], y[5])
        rhsx11 = diffx(A, x[5], y[5])
        rhsy12 = diffy(A, x[5], y[5])

        phis7 = phis(x[3], y[3])
        phisx8 = phisdiffx(phisbasis, x[3], y[3])
        phisxx9 = phisdiffxx(phisbasis, x[3], y[3])
        phisy9 = phisdiffy(phisbasis, x[3], y[3])
        phis10 = phis(x[5], y[5])
        phisx11 = phisdiffx(phisbasis, x[5], y[5])
        phisy12 = phisdiffy(phisbasis, x[5], y[5])

        Amat[7,:] = phis7
        Amat[8,:] = phisx8
        Amat[9,:] = phisxx9 + N3*phisy9
        Amat[10,:] = phis10
        Amat[11,:] = phisx11
        Amat[12,:] = phisy12

        b = [rhs1, rhsy2, rhsyy3+N1*rhsx3, rhs4, rhsy5, rhsyy6+N2*rhsx6,
        rhs7, rhsx8, rhsxx9+N3*rhsy9, rhs10, rhsx11, rhsy12]
    else
        rhs7 = calcphi(A, x[4], y[4])
        rhsx8 = diffx(A, x[4], y[4])
        rhsy9 = diffy(A, x[4], y[4])
        rhs10 = calcphi(A, x[5], y[5])
        rhsx11 = diffx(A, x[5], y[5])
        rhsy12 = diffy(A, x[5], y[5])

        phis7 = phis(x[4], y[4])
        phisx8 = phisdiffx(phisbasis, x[4], y[4])
        phisy9 = phisdiffy(phisbasis, x[4], y[4])
        phis10 = phis(x[5], y[5])
        phisx11 = phisdiffx(phisbasis, x[5], y[5])
        phisy12 = phisdiffy(phisbasis, x[5], y[5])

        Amat[7,:] = phis7
        Amat[8,:] = phisx8
        Amat[9,:] = phisy9
        Amat[10,:] = phis10
        Amat[11,:] = phisx11
        Amat[12,:] = phisy12

        b = [rhs1, rhsy2, rhsyy3+N1*rhsx3, rhs4, rhsy5, rhsyy6+N2*rhsx6,
        rhs7, rhsx8, rhsy9, rhs10, rhsx11, rhsy12]
    end

    c = Amat \ -b
    phi = phifuncts()
    particularphi = calcphifunct(A)
    f(x,y) = particularphi(x,y) .+ dot(c, [p(x,y) for p in phi])

    return f
end

psi = gen_flux_function(A, delta, epsilon, kappa, [1-1.1*delta*epsilon, 1.1*kappa*epsilon], [0,0])
for i=1:200
    for j=1:200
        z[j,i] = psi(x[i], y[j])
    end
end

contour(z, levels=100)
