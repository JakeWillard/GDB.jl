
import numpy as np
import sympy as sp


def cerfon_flux(x, y, A, C):
    """
    Defines general solution to Grad-Shafranov from Cerfon et al 2010.
    Defines a function of x and y with parameters A and C according to the
    expansion detailed in Cerfon's paper which solves the Grad-Shafranov equation
    with Solov'ev profiles.
    Args:
        x (symbol): x coordinate
        y (symbol): y coordinate
        A (symbol): constant that determines the B_phi profile
        C (indexed symbol): coefficients for expansion
    Returns:
        expression: general solution for flux function
    """

    # define vector for homogenous basis function
    psi_H = sp.zeros(1, 12)

    # define particular and homogenous terms
    psi_p = (1 - A) * x**4 / 8 + A * x**2 * sp.log(x) / 2
    psi_H[0, 0] = 1
    psi_H[0, 1] = x**2
    psi_H[0, 2] = y**2 - x**2 * sp.log(x)
    psi_H[0, 3] = x**4 - 4 * x**2 * y**2
    psi_H[0, 4] = 2 * y**4 - 9 * y**2 * x**2 + 3 * x**4 * \
        sp.log(x) - 12 * x**2 * y**2 * sp.log(x)
    psi_H[0, 5] = x**6 - 12 * x**4 * y**2 + 8 * x**2 * y**4
    psi_H[0, 6] = 8 * y**6 - 140 * y**4 * x**2 + 75 * y**2 * x**4 - 15 * x**6 * \
        sp.log(x) + 180 * x**4 * y**2 * sp.log(x) - \
        120 * x**2 * y**4 * sp.log(x)
    psi_H[0, 7] = y
    psi_H[0, 8] = y * x**2
    psi_H[0, 9] = y**3 - 3 * y * x**2 * sp.log(x)
    psi_H[0, 10] = 3 * y * x**4 - 4 * y**3 * x**2
    psi_H[0, 11] = 8 * y**5 - 45 * y * x**4 - 80 * y**3 * \
        x**2 * sp.log(x) + 60 * y * x**4 * sp.log(x)

    # return general solution
    return psi_p + np.sum([C[i] * psi_H[0, i] for i in range(0, 12)])


def fit_geometry(x, y, psi, C, ep, de, ka, r, x_up, x_down):
    """
    Places geometric constraints on psi and solves for coefficients C.
    Fixes value and derivatives of psi at specific points to determine
    coefficients C that produce a desired flux function. This essentially is
    doing the same calculations discussed in Cerfon et al 2010. Returns a list of
    sympy substutution tuples for expressions that depend on C.
    Args:
        x (symbol): x coordinate
        y (symbol): y coordinate
        psi (expression): general solution for flux function, depends on (x, y, C)
        C (indexed symbol): coefficients for general solution
        ep (float): aspect ratio
        de (float): triangulation
        ka (float): elongation
        r (float): shift parameter that locates the x-points
        x_up (boolean): True for upper x-point, False for no upper x-point
        x_down (boolean): True for lower x-point, False for no lower x-point
    Returns:
        list: tuples (C[i], C_solved[i]) for each i
    """

    # compute x-point shifts
    dx = de * ep * (1 + r)
    dy = ka * ep * (1 + r)

    # compute subs lists for points of interest
    a = [(x, 1 + ep), (y, 0)]
    b = [(x, 1 - de * ep), (y, ka * ep)]
    c = [(x, 1 - ep), (y, 0)]
    d = [(x, 1 - de * ep), (y, -ka * ep)]
    xp1 = [(x, 1 - dx), (y, dy)]
    xp2 = [(x, 1 - dx), (y, -dy)]

    # compute N values
    alpha = np.arcsin(de)
    N1 = -(1 + alpha)**2 / (ep * ka**2)
    N2 = (1 - alpha)**2 / (ep * ka**2)
    N3 = -ka / (ep * np.cos(alpha)**2)

    # initialize list of constraints
    eqns = []

    # add inboard and outboard side constraints
    eqns.append(psi.subs(a))
    eqns.append(psi.diff(y, 1).subs(a))
    eqns.append((psi.diff(y, 2) + N1 * psi.diff(x, 1)).subs(a))
    eqns.append(psi.subs(c))
    eqns.append(psi.diff(y, 1).subs(c))
    eqns.append((psi.diff(y, 2) + N2 * psi.diff(x, 1)).subs(c))

    # If upper x-point, set B=0 at x-point. Otherwise,
    # set direction and curvature of core-side surface
    if x_up:
        eqns.append(psi.subs(xp1))
        eqns.append(psi.diff(x, 1).subs(xp1))
        eqns.append(psi.diff(y, 1).subs(xp1))
    else:
        eqns.append(psi.subs(b))
        eqns.append(psi.diff(x, 1).subs(b))
        eqns.append((psi.diff(x, 2) + N3 * psi.diff(y, 1)).subs(b))

    # If lower x-point, set B=0 at x-point. Otherwise,
    # set direction and curvature of core-side surface
    if x_down:
        eqns.append(psi.subs(xp2))
        eqns.append(psi.diff(x, 1).subs(xp2))
        eqns.append(psi.diff(y, 1).subs(xp2))
    else:
        eqns.append(psi.subs(d))
        eqns.append(psi.diff(x, 1).subs(d))
        eqns.append((psi.diff(x, 2) - N3 * psi.diff(y, 1)).subs(d))

    # solve for C
    Cs = list(sp.linsolve(eqns, [C[i] for i in range(0, 12)]))[0]

    # return list of substitutions for C
    return [(C[i], Cs[i]) for i in range(0, 12)]
