# ImmersedMirrors.jl Guide

### Dependencies

* ProgressMeter.jl
* ForwardDiff.jl
* RecipesBase.jl

### Overview

This is a module for setting up immersed boundary problems in 2D using three basic structs: Mirror, Grid, and GhostData.

A Mirror struct is used to represent a polygonal boundary for multiple-connected domains. As the name suggests, this boundary will be treated as a kind of mirror surface, so that points can be reflected across the boundary. The reflection operations are key to how we handle immersed boundary conditions.

A Grid struct represents the points used to discretize the space. Grid-points are arranged on a square lattice, but the overall arrangement is in general non-rectangular. Rather, the grid-points are chosen not to cover much more than the region enclosed by the boundaries, with a variable excess of exterior "ghost-points" which are needed to properly enforce boundary conditions.

A GhostData struct contains the output of an important preprocessing step in our approach to immersed boundary problems. It is with this struct that the values of quantities on ghost-points can be extrapolated from interior values in accordance with symmetric or anti-symmetric conditions. This functionality is all that is needed to enforce Neumann-like and Dirichlet-like boundary conditions.

### Mirror

#### Constructor

A Mirror struct is initialized with a single input having type Vector{Matrix}. Each element of the vector is a 2xN matrix, N being some number, representing a closed chain of points (the first vertex of each chain is not double-counted). As the name suggests, these chains are meant to represent mirror surfaces in the plane: oriented curves that points in space can be reflected about. The following code creates a mirror out of 2 different circles approximated by 10-sided regular polygons:

    thetas = LinRange(0, 2*pi, 11)[1:10]
    chain1 = zeros(2, 10)
    chain2 = zeros(2, 10)
    for i=1:10
        chain1[:,i] = [cos(thetas[i]) - 1, sin(thetas[i])]
        chain2[:,i] = [cos(thetas[i]) + 1, sin(thetas[i])]
    end

    M = Mirror([chain1, chain2])

The orientation of the chains will be important, so note how the vertices are listed in counter-clockwise order for both circles. Since all chains are closed loops, the orientation determines which side of the chain is considered "inside" or "outside". For vertices listed in counter-clockwise order, points enclosed by the loop are "inside" points. If we reverse the order, that will also reverse the orientation, so the points enclosed by the loop will be "outside" points. We choose the orientation based on whether we want the chain the represent an outer boundary of a domain, or the boundary of a "hole" in the domain. For example, we can define a Mirror corresponding with an annalus:

    thetas = LinRange(0, 2*pi, 11)[1:10]
    chain1 = zeros(2, 10)
    chain2 = zeros(2, 10)
    for i=1:10
        chain1[:,i] = 2*[cos(thetas[i]), sin(thetas[i])]
        chain2[:,i] = [cos(thetas[i]), -sin(thetas[i])]
    end

    M = Mirror([chain1, chain2])

#### Usage

distance_to_mirror(x, y, M::Mirror) returns a tuple (dist, k) where dist is the shortest distance from the point (x,y) to the mirror struct M, and k is the index corresponding to the nearest vertex. All of the vertices of all of the chains are listed together in the attribute "verts", so this nearest vertex can be referenced with M.verts[:,k]. The returned value dist is a signed distance, and it returns positive for interior points and negative for exterior points.

mirror_image(x, y, M::Mirror) maps exterior points (x,y) to interior points with respect to a Mirror struct M, in a way corresponding to a reflection across the nearest chain. The function acts as an identity operation for interior points, i.e. (x,y) = mirror_image(x, y, M) if distance_to_mirror(x, y, M)[1] > 0.

### Grids

#### Constructor

A Grid struct is initialized with three pieces of info: (1) a Mirror struct to define the boundary for our problem, (2) a parameter determining the number of ghost-points, (3) a bounding box that contains the domain, and (4) a parameter determining the resolution of the grid. The function signature is therefore Grid(M, h, r0, r1, p), where the input variables are:
* M -- a Mirror struct
* h -- thickness of the ghost region. (A ghost point satisfies -h < distance_to_mirror(x,y,M) < 0)
* r0 -- vector for the bottom left corner of the bounding box.
* r1 -- vector for the top right corner of the bounding box.
* p -- N=2^p, where N is the number of grid-points along the shortest side of the bounding box. For example, if (r1 - r0)[1] < (r1 - r0)[2], then Nx = 2^p. If (r1 - r0)[1] > (r1 - r0)[2], then Ny = 2^p.

#### The points attribute and function_to_grid

Grids have an attribute called "points", which is of size 2 x Nk, which Nk is the number of grid-points. Each column of points is a position (x,y) inside the domain. We can use Grids to evaluate a function at each of these points and store those variables in a vector. We do this with function_to_grid:

    f_as_vector = function_to_grid(grd::Grid) do x,y
        sin(x)*sin(y)
    end

    f_as_vector_alt = []
    for k=1:grd.Nk
        x, y = grd.points[:,k]
        f_as_vector = [f_as_vector_alt; sin(x)*sin(y)]
    end

    @assert all(f_as_vector .== f_as_vector_alt)

The above code demonstrates the usage of function_to_grid, while also showing what that function does by constructing the output vector with an explicit loop over grd.points.

#### The projection matrix

It is convenient to think of vectors like f_as_vector defined above as living in a vector space associated with the list of points grd.points. So from here on, I will talk about a "grid space" as being that kind of vector space. Remember that the Grid constructor used input data corresponding to a rectangular grid, with corners at r0 and r1, which has some dimension Nx x Ny. In general, this ends up being translated into a larger rectangular grid, with dimension _Nx x _Ny where _Nx = Nx + 2Nbuffer and _Ny = Ny + 2Nbuffer, Nbuffer is an optional keyword argument with default value 100. This larger rectangle also has different corners offset by Nbuffer*grd.dr. This rectangular grid also has a vector space, and it is convenient to be able to define certain things on the rectangular space and then project them onto the non-rectangular space. The Proj attribute is a sparse matrix that does exactly this. The following should make it clear:

    fvec_non_rec = function_to_grid(grd::Grid) do x,y
        sin(x)*sin(y)
    end

    fvec_rec = Float64[]
    for j=1:grd._Ny
        for i=1:grd._Nx
            x = grd.r0[1] + grd.dr*(i-1)
            y = grd.r0[2] + grd.dr*(j-1)
            fvec_rec = [fvec_rec; sin(x)*sin(y)]
        end
    end

    @assert all(fvec_non_rec .== grd.Proj * fvec_rec)

Here, rectangular and non-rectangular version of the function sin(x)sin(y) are calculated, and it is asserted that the non-rectangular version can be obtained from the rectangular version by applying the matrix transformation grd.Proj.

#### Mapping linear operators

It is often convenient to define linear operators on the rectangular space and then "translate" that matrix onto the non-rectangular space. This can only be done exactly if the transformation is invertible, and Proj is not generally invertible. However, there is an appropriate pseudoinverse in this case, which is Proj^T. This can be seen as a mapping from the non-rectangular space to the rectangular space assuming that everything is zero at points not contained in grd.points. This is an acceptable assumption to make, provided that the "ghost region" is sufficiently thick. In other words, h/grd.dr must not be too small, where h is the input variable controlling the thickness of the ghost region. If h and the grid resolution are chosen carefully, linear operators defined on the rectangular grid can be mapped onto the non-rectangular grid via Proj M Proj^T, and the resulting matrix will be a faithful representation of the matrix M everywhere inside the domain of interest.

operator_to_grid is a function that will automatically project a matrix M defined on the rectangular space onto the non-rectangular space using Proj M Proj^T. This is done by defining M within a do-block:

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

The matrix Diff_x is just a first derivative matrix calculated using a central difference approximation. Note that the code within the do-block is nothing more than a crude method for computing a central difference matrix on a rectangular grid.


### GhostData

#### Constructor

GhostData structs are a marriage between a Mirror and Grid struct. They use the data composed in each to carry out a preprocessing step useful in solving boundary-value problems. It is advisable to run the constructor using multiple processors since it can involve a lot of calculation. The constructor can be called simply with GhostData(M::Mirror, grd::Grid).

#### Ghost value extrapolation

To impose boundary conditions, we need to set values at ghost points according to values at interior points. For Neumann boundary conditions, we can do this with a symmetric condition: apply the mirror_image function to a ghost point, and from that mirror image find the nearest grid-point enclosed by the boundary, then the symmetric condition is applied by setting the ghost value equal to the mirror value. For Dirichlet conditions, we instead apply anti-symmetric conditions, where the ghost value is set equal to minus the mirror value. Both of these operations are linear, and therefore can be represented by a matrix transformation. The GhostData constructor creates this matrix as an attribute called "R".

By default, R will always apply a symmetric extrapolation. To get anti-symmetric extrapolation, we use the flip_segment function. Given a list of indices inds corresponding to specific segments of the polygon boundary, we can "flip" between symmetric and anti-symmetric conditions with respect to those specific segments. flip_segment(gd::GhostData, inds) returns a new GhostData struct with a modified R matrix, where anti-symmetric conditions will now be used for the segments corresponding with inds. flip_segment will also turn anti-symmetric conditions into symmetric conditions, so applying the function again to the new GhostData struct will result in a copy of the original GhostData struct we started with.

#### Imposing boundary conditions (simple)

Say we want to numerically evolve the Schroedinger equation for a free particle. Having defined a GhostData struct called gd, we could easily account for homogeneous Dirichlet or Neumann conditions like this:

    # define laplacian matrix
    L = operator_to_grid(grd::Grid) do
        Dxx = spdiagm(1 => ones(grd._Nx-1))
        Dxx += spdiagm(-1 => ones(grd._Nx-1))
        Dxx += spdiagm(0 => -2*ones(grd._Nx))

        Dyy = spdiagm(1 => ones(grd._Ny-1))
        Dyy += spdiagm(-1 => ones(grd._Ny-1))
        Dyy += spdiagm(0 => -2*ones(grd._Ny))

        kron(Dyy, Dxx) / grd.dr^2
    end

    # update wavefunction with forward-Euler (hbar/m = 2)
    psi += dt*im*L*psi

    # impose boundary conditions
    psi = gd.R*psi

The code following the definition of the laplacian would be all it takes to update the wavefunction by one timestep while satisfying the boundary conditions, provided that forward-Euler method was deemed appropriate.

We can modify this to satisfy inhomogeneous conditions as well. This is done by instead imposing (psi - psi_b) = gd.R*(psi - psi_b), where psi_b is either the boundary value for a Dirichlet condition, or its normal derivative at the boundary is the boundary value for a Neumann condition.

#### Constraining linear systems

Attached to GhostData structs is an attribute called Proj, which is another projection matrix that is similar but distinct from the Proj attached to Grid structs. Here, Proj is a projection from the non-rectangular space to an even smaller subspace defined by the points in grd.points that are not ghost points. From here, I will call the non-rectangular space that includes grid-points as the ghost-inclusive space, and this smaller subspace the ghost-exclusive space.

The purpose of defining yet another projection onto yet another grid-space is clear if we make the following calculation:

    # function returns laplacian on the ghost-inclusive space
    L = operator_to_grid(grd::Grid) do

        # code in the block creates laplacian on rectangular space
        Dxx = spdiagm(1 => ones(grd._Nx-1))
        Dxx += spdiagm(-1 => ones(grd._Nx-1))
        Dxx += spdiagm(0 => -2*ones(grd._Nx))

        Dyy = spdiagm(1 => ones(grd._Ny-1))
        Dyy += spdiagm(-1 => ones(grd._Ny-1))
        Dyy += spdiagm(0 => -2*ones(grd._Ny))

        kron(Dyy, Dxx) / grd.dr^2
    end

    # one more projection to get laplacian on ghost-exclusive space.
    L = gd.Proj * L * gd.R * transpose(gd.Proj)

This code calculates a laplacian matrix again, but then carries out another projection of the form Proj L Proj^(-1). However, instead of using Proj^T as the pseudoinverse, now we use R Proj^T. R Proj^T will map vectors on the ghost-exclusive space into the ghost-inclusive space such that the ghost values are extrapolated via the transformation described previously. In this way, L is now a matrix representation of the laplacian that acts on the ghost-exclusive space in a way that automatically imposes homogeneous boundary conditions. This means that we can solve simple PDEs:

    rhs = function_to_grid(grd) do x,y
        sin(x)*sin(y)
    end
    rhs = gd.Proj * rhs

    # solve poisson equation on ghost-exclusive space
    phi_projected = L \ rhs

    # extrapolate with gd.R to get solution on ghost-inclusive space
    phi = gd.R * transpose(gd.Proj) * phi_proj

We have just solved the Poisson equation Del^2 phi = sin(x)sin(y) with arbitrary homogeneous boundary conditions defined by the GhostData struct gd.

#### Inhomogeneous generalization

The fully generalized version of the above procedure for inhomogeneous boundary conditions is wrapped up in the function constrain_system(A::SparseMatrixCSC, b::Vector{Float64}, xb::Vector{Float64}, gd::GhostData). A and b represent the linear problem we want to solve Ax=b, and xb is the vector representing the function that defines our boundary value.

constrain_system returns a tuple (Anew, bnew, Pi, c), where we obtain our solution on the ghost-exclusive space by solving Anew x = bnew, and then that is mapped onto the ghost-inclusive space via Pi x + c. The following repeats the calculation made in the previous section, but with an inhomogeneous boundary condition:

    rhs = function_to_grid(grd) do x,y
        sin(x)*sin(y)
    end

    fb = function_to_grid(grd) do x,y
        norm([x,y])
    end

    Anew, bnew, Pi, c = constrain_system(L, rhs, fb, gd::GhostData)
    phi = Pi * (Anew \ bnew) + c
