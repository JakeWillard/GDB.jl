

## ImmersedMirrors.jl Quick Guide

### Dependencies

* ProgressMeter.jl
* ForwardDiff.jl
* RecipesBase.jl


### Mirrors

#### Define a Mirror

All you need for this is to represent the boundary of the domain as a collection of polygons. Each polygon should be represented by a list of vertices, and the order in which they are listed determines what regions get excluded from the domain: counter-clockwise results in the exterior being excluded from the domain, and clockwise results in the interior being excluded. This allows for the possibility of multiply-connected domains.

Call the constructor with Mirror(chains::Vector{Matrix{Float64}}). The input variable chains is a list of matrices, each matrix having size 2xN for some N, where N is the number of vertices for the particular chain. Here is an example:

    thetas = LinRange(0, 2*pi, 11)[1:10]
    chain1 = zeros(2, 10)
    chain2 = zeros(2, 10)
    for i=1:10
      chain1[:,i] = 2*[cos(thetas[i]), sin(thetas[i])]
      chain2[:,i] = [cos(thetas[i]), -sin(thetas[i])]
    end

    M = Mirror([chain1, chain2])

This code creates a Mirror struct for a domain approximating an annalus with inner radius 1 and outer radius 2.


#### Distances and mirror images


You can compute a signed distance to the boundary using distance_to_mirror(x,y,M::Mirror), which returns a tuple (dist, k) where dist is the signed distance from the point (x,y) to the boundary represented by the Mirror struct M, and k is an index for the nearest vertex, which you can subsequently reference via M.verts[:,k]. The distance is signed because it is positive or negative depending on the orientation of the polygon. The distance is negative for points excluded from the domain, and it is positive for points included in the domain.

As the name suggests, a Mirror struct is meant to represent the boundary as a reflecting surface. To reflect a point (x,y) across the boundary, use mirror_image(x,y,M::Mirror). However, a Mirror in this case is only reflecting for points excluded from the domain, since we use this to reflect ghost-points. For included points the function acts as an identity, i.e. (x,y) = mirror_image(x, y, M) if distance_to_mirror(x, y, M)[1] > 0.


### Grids

#### Defining a Grid

The Grid struct contains information associated with a non-rectangular grid with a square lattice that covers the domain defined by a Mirror struct. It is non-rectangular in general, since the grid-points will be arranged in a way to cover not much more than the points included in the domain. The grid will include ghost-points outside the domain within some given distance, so that the boundary is immersed in the grid and has some number of points surrounding it.

Use the constructor Grid(M::Mirror, h::Float64, r0::Vector{Float64}, r1::Vector{Float64}, p::Int64). The meanings of the inputs are:
* M is the mirror struct that defines the domain.
* h is the thickness of the ghost region.
* r0 is the bottom-left position of the bounding box.
* r1 is the top-right position of the bounding box.
* p determines the resolution: # of grid-points ~ 2^(2p).

r0 and r1 determine a bounding box that begins the arrangement of the grid-points. The arrangement starts off as a rectangular grid with dimension Nx x Ny where r0 and r1 are the bottom-left and top-right corners respectively. p determines the resolution by setting the number of grid-points along the shortest dimension: if (r1 - r0)[2] > (r1 - r0)[1], then Nx = 2^p, but if (r1 - r0)[2] < (r1 - r0)[1], then Ny = 2^p. The lattice is square, so this calculation determines the spacing between grid-points dr. The number of grid-points in the other direction is then fixed to be whatever fills the bounding box.

The following code creates a grid to cover a domain defined by a Mirror struct M:

    r0 = [-1.0, -2.0]
    r1 = [1.0, 2.0]
    h = 0.1
    grd = Grid(M::Mirror, h, r0, r1, 8)

The bounding box for this region is a 2 x 4 rectangle, so Nx will be 2^8, the spacing will be dr=2/Nx, and Ny will be div(4, dr). To see what the actual grid-points are, you can reference them with grd.points, which is a 2xNk matrix where Nk is the number of grid-points (always less than or equal to Nx*Ny). Nk can be referenced as well with grd.Nk.

#### Representing and plotting functions

Functions of (x,y) have a representation on the grid as vectors of length Nk. You can easily compute such vectors:

    f(x,y) = sin(x)*sin(y)
    f_as_vector = function_to_grid(f, grd::Grid)

    @assert all(f_as_vector .== [f(grd.points[:,k]...) for k=1:grd.Nk])

This code created a vector representing the function sin(x)sin(y). The last part of the code shows what function_to_grid does: it is nothing more than the calculation [f(grd.points[:,k]) for k=1:grd.Nk].

This module includes a plotting recipe for vectors representing functions, which will make Plots-compatible function interpret the vector as 2D data:

    f_as_vector = function_to_grid(grd::Grid) do x,y
        sin(x)*sin(y)
    end

    using Plots
    heatmap(f_as_vector, grd::Grid)
    contour(f_as_vector, grd::Grid)

I repeated the calculation of f_as_vector here to illustrate the do-block syntax.


#### Representing linear operators

Linear operators that act on functions of (x,y) also have a representation of the grid as matrices. Computing matrices like this might be tricky to do since the grid is non-rectangular in general. For things like finite difference approximations, it would be much easier if the grid were rectangular at least. The Grid struct offers a convenient way to compute matrix representations for the grid provided that you know how to compute them on a non-rectangular grid:

    Diff_x = operator_to_grid(grd::Grid) do
        Dx = spdiagm(1 => ones(grd._Nx-1))
        Dx += spdiagm(-1 => -ones(grd._Nx-1))
        Iy = sparse(I, grd._Ny, grd._Ny)
        kron(Iy, Dx) / grd.dr
    end

    fvec = function_to_grid(grd::Grid) do x,y
        sin(x)*sin(y)
    end

    fvec_diffx = Diff_x * fvec

Diff_x is a matrix representing a partial derivative with respect to x using a central different approximation. This matrix is used to compute fvec_diffx, which is a vector approximating the partial derivative of sin(x)sin(y) with respect to x. Note that the code inside the do-block is nothing more than a crude method for constructing a central difference matrix on a rectangular grid with dimension _Nx x _Ny. _Nx and _Ny are attributes of the Grid struct, and are not exactly the same as Nx and Ny described before. What they exactly correspond to is technical, the only thing we need to know here is that the code in the do-block must assume a rectangular grid having grd._Nx and grd._Ny as dimensions, and a square lattice with spacing grd.dr. The rest is taken care of by the function operator_to_grid, which maps the rectangular version of the operator to a suitable matrix representation for our grid.

One important note is that this calculation relies on finding a mapping between the grid and a rectangular grid where the code in the do-block applies, which isn't actually possible to do exactly without some implicit assumption. Here, the implicit assumption will always be that functions are zero in fully excluded regions of the x-y plane, i.e. points that have distance_to_mirror(x,y,M)[1] < -h. This could be unacceptable in some cases, but for matrices that merely do local operations (finite difference matrices) the error will be localized to the very edge of the grid. So in those cases, you just need to be aware that the matrices won't be accurate representation of the operator on all of the ghost points. This won't be a problem in the case of boundary value problems using the extrapolation method described in the next section, and provided that the ghost region is not made too small.


### Extrapolators


#### Defining an Extrapolator

An Extrapolator struct just stores the output of an important pre-processing step when setting up boundary value problems. All you need to do the calculation is call Extrapolator(M::Mirror, grd::Grid). This calculation could take some time, especially at high resolutions, so it is advisable to run this with multiple processors (does not necessitate shared memory).


#### Using the Extrapolator for Neumann conditions.

By default, an Extrapolator can be used to impose Neumann boundary conditions (not necessarily homogeneous) by extrapolating values at ghost-points according to a symmetry condition. This is where the mirror_image function comes in: for homogeneous Neumann conditions, the value at a ghost-point (x,y) is simply set equal to the interior value found at the ghost point's mirror image using mirror_image(x,y,M::Mirror).

To transform a vector to alter the ghost values in this way, simply pass the vector into the Extrapolator struct itself:

    fvec = function_to_grid(grd::Grid) do x,y
        sin(x)*sin(y)
    end

    extr = Extrapolator(M::Mirror, grd::Grid)

    boundary_value = zeros(grd.Nk)

    fvec_ex = extr(fvec, boundary_value)

Note that the Extrapolator instance is callable like a function. Here, fvec_ex is the same as fvec, but components corresponding to ghost points have be set according to the symmetry condition, which will set the derivative at the boundary equal to the derivative of the function represented by boundary_value. For this example, boundary_value is set to zero, corresponding to a homogeneous boundary condition.

A linear system can be constrained to implicitly satisfy the boundary conditions by passing different info into the Extrapolator:

    # compute extrapolator
    extr = Extrapolator(M::Mirror, grd::Grid)

    # define laplacian matrix
    L = operator_to_grid(grd) do
        Dxx = spdiagm(1 => ones(grd._Nx-1))
        Dxx += spdiagm(-1 => ones(grd._Nx-1))
        Dxx += spdiagm(0 => -2*ones(grd._Nx))

        Dyy = spdiagm(1 => ones(grd._Ny-1))
        Dyy += spdiagm(-1 => ones(grd._Ny-1))
        Dyy += spdiagm(0 => -2*ones(grd._Ny))

        Ix = sparse(I, grd._Nx, grd._Nx)
        Iy = sparse(I, grd._Ny, grd._Ny)

        (kron(Iy, Dxx) + kron(Dyy, Ix)) / grd.dr^2
    end

    # right hand side vector
    rhs = zeros(grd.Nk)

    # boundary value
    bdry = function_to_grid(grd::Grid) do x,y
        sin(x)^2
    end

    A, b = extr(L, rhs, bdry)
    phi = extr(A \ b, bdry)

This code solves Laplace equation Del^2 phi = 0 with inhomogeneous Neumann condition, where the normal derivative of phi - sin^2(x) is zero on the boundary. Calling the Extrapolator extr with the signature extr(L::SparseMatrixCSC, rhs::Vector{Float64}, bdry::Vector{Float64}) will return a new sparse matrix A and vector b. The linear system A x = b is a lower dimensional version of the system L x = rhs, where the boundary conditions have been automatically imposed on the system. Because the system is lower dimension, one would rightly expect that solving it via A \ b does not immediately result in a vector that represents the solution on our grid. To get that, we pass it into the Extrapolator again along with the boundary value.

#### Dirichlet conditions

You can also generate Extrapolators that impose anti-symmetry conditions instead of symmetry conditions. The calculations involved are very similar to how the symmetry conditions are imposed, so there is no need to redo the preprocessing calculations. The change is made on a by segment basis, i.e. it is possible to impose a mix of Dirichlet and Neumann conditions on different parts of the boundary, and each side of the boundary polygon can be independently switched to Dirichlet conditions.

    extr0 = Extrapolator(M, grd)

    inds = 1:5
    extr = make_dirichlet(extr0, inds)

In this case, we took an Extrapolator extr0 and switched the first 5 segments of the boundary (corresponding to the 5 vertices M.verts[:,1:5]) so that extr will now extrapolate ghost points near those segments using anti-symmetry conditions. Its not clear what this corresponds with geometrically, since we never defined the surface corresponding to the Mirror struct. Here is a more clear example:

    # make Mirror corresponding to annalus
    thetas = LinRange(0, 2*pi, 11)[1:10]
    chain1 = zeros(2, 10)
    chain2 = zeros(2, 10)
    for i=1:10
        chain1[:,i] = 2*[cos(thetas[i]), sin(thetas[i])]
        chain2[:,i] = [cos(thetas[i]), -sin(thetas[i])]
    end
    M = Mirror([chain1, chain2])

    # define Grid
    r0 = [-2.1, -2.1]
    r1 = [2.1, 2.1]
    h = 0.1
    grd = Grid(M::Mirror, h, r0, r1, 8)

    # default Neumann Extrapolator
    extr_neumann = Extrapolator(M, grd)

    # make outer boundary dirichlet
    inds = 1:10
    extr = make_dirichlet(extr_neumann, inds)

Here, we defined the boundary to approximate an annalus like we did at the beginning of this guide. The vertices M.verts[:,1:10] correspond to the outer boundary, so we can make the Extrapolator impose Dirichlet conditions on the outer boundary by passing 1:10 into the make_dirichlet function.

Even though it might seem as though there are two distinct surfaces here (the inner and outer boundary), it's important to note that the boundary conditions are not imposed seperetly. We still define a single vector for the boundary value, and extrapolation operations are evaluated together. As an example, here is how we could impose a homogeneous Neumann condition at the inner boundary and an inhomogeneous Dirichlet condition at the outer boundary:

    bdry = function_to_grid(grd::Grid) do x,y
        r = norm([x,y])
        if r < 1.5
            return 0.0
        else
            return r
        end
    end

    A, b = extr(L, zeros(grd.Nk), bdry)
    phi = extr(A \ b, bdry)

This code solves the Laplace equation again, this time with a boundary value vector that is zero everywhere that r < 1.5, which will result in a homogeneous Neumann condition at the inner boundary, and equal to r for r > 1.5, which will result in a Dirichlet condition phi = r at the outer boundary.
