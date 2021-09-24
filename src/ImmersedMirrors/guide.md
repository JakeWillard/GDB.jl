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

For brevity, I will from here on refer to counter-clockwise ordering as "positive" orientation and clockwise as "negative" orientation. Positive orientation defines bounded regions, negative orientation defines unbounded regions with a hole.

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

### The points attribute and function_to_grid

Grids have an attribute called "points", which is of size 2 x Nk, which Nk is the number of grid-points. Each column of points is a position (x,y) inside the domain. We can use Grids to evaluate a function at each of these points and store those variables in a vector. We do this with function_to_grid:

    f_as_vector = function_to_grid(grd::Grid) do x,y
        sin(x)*sin(y)
    end

    f_as_vector_alt = []
    for k=1:grd.Nk
        x, y = grd.points[:,k]
        f_as_vector = [f_as_vector; sin(x)*sin(y)]
    end

    @assert all(f_as_vector .== f_as_vector_alt)

The above code demonstrates the usage of function_to_grid, while also showing precicely what that function does by constructing the output vector with an explicit loop over grd.points.

### The projection matrix

Remember that the Grid constructor used input data corresponding to a rectangular grid, with corners at r0 and r1, which has some dimension Nx x Ny. Having a vector of data like f_to_vector computed in the above code, it is sometimes useful to have this in an alternate representation corresponding to the Nx x Ny rectangular grid instead. It is easiest to explain this by explicit example:

    grd = Grid([-1.0, -1.0], [1.0, 1.0], 0.01, 1) do x,y
        norm([x,y]) < 1
    end

    fvec_non_rec = function_to_grid(grd::Grid) do x,y
        sin(x)*sin(y)
    end

    fvec_rec = []
    for j=1:grd._Ny
        for i=1:grd._Nx
            x = -1.0 + 0.01*(i-1)
            y = -1.0 + 0.01*(j-1)

            if norm([x,y]) < 1
                fvec_rec = [fvec_rec; sin(x)*(sin(y))]
            else
                fvec_rec = [fvec_rec; 0]
            end
        end
    end

    @assert all(fvec_non_rec .== grd.Proj*fvec_rec)
    @assert all(fvec_rec .== transpose(grd.Proj)*fvec_non_rec)

The first part of the code constructs a Grid. The bounding box is a square centered at the origin with side length 2, the lattice length 0.01, and the domain is defined as the unit disk. The second part of the code uses function_to_grid to calculate sin(x)sin(y) at each of the points inside the unit disk. The third part of the code makes a similar calculation, but instead of only storing values at the points contained in the grd.points attribute, each point on the Nx x Ny grid has a corresponding entry in fvec_rec, and points not within the unit disk are assigned a value of zero. These two vectors fvec_non_rec and fvec_rec are related by linear transformation given by the attribute grd.Proj, which is a sparse matrix. The final part of the code asserts this identity, which is that grd.Proj transforms the rectangular representation into the non-rectangular representation, and its transpose is a right-inverse.

### operator_to_grid

Linear operators should have some representation on the grid as well. The derivative of sin(x)sin(y), for example, should be approximated by some matrix transformation of f_as_vector defined above. In the case of a rectangular grid, such matrices are very easy to compute. The utility of the projection matrix described in the previous section is that it allows us to compute linear operators for the non-rectangular case using rectangular versions. The vector representation of the non-rectangular grid is basically a subspace of the rectangular grid representation, and so the transformation of linear operators works exactly the same as a change of basis. Given a matrix transformation computed for the rectangular grid (size NxNy x NxNy), the non-rectangular version is simply Proj M Proj^T. The function operator_to_grid helps streamline the calculation of matrix operators:

    Diff_x = operator_to_grid(grd) do
        Dx = spdiagm(1 => ones(grd._Nx-1))
        Dx += spdiagm(-1 => -ones(grd._Nx-1))
        Iy = sparse(I, grd._Ny, grd._Ny)
        kron(Iy, Dx) / grd.dr
    end

    fvec = function_to_grid(grd::Grid) do x,y
        sin(x)*sin(y)
    end

    fvec_diffx = Diff_x * fvec

The code within the do-block is nothing more than a crude method of calculating a central difference approximation of the first derivative assuming a rectangular grid. It is automatically projected onto the non-rectangular subspace, so that fvec_diffx is an approximation of the partial x derivative of sin(x)sin(y) on the points contained in grd.points.

### Projection errors

The matrix Diff_x computed above has some potential errors associated with it. The first is due to the crudeness of implementation: the code within the do-block is simply not a correct central difference approximation at the edges of the rectangular grid. This is purely user-error, since I could have written every row of the matrix to be computed correctly. However, we can get away with sloppiness of this kind provided that the edges of the rectangular grid are not very near the edges of the actual domain, i.e the bounding box defined by r0 and r1 should be large enough so that these errors are "deleted" by the similarity transform Proj M Proj^T.

The other source of error has to do with the assumptions implicit in the projection itself. Proj is not invertible, rather we are using Proj^T as the appropriate pseudoinverse: Proj^T b is the minimum norm solution of Proj x = b. In the case of operators which couple points inside the domain to points outside the domain, Proj M Proj^T will not give the same result as M in general, as it effectively will perform the transformation under the assumption that everything outside the domain is set to zero. Tools described later in this guide will be useful in properly accounting for this effect when solving boundary value problems.


### GhostData

#### Constructor

A GhostData struct exists as a marriage between a Mirror and a Grid, and it is initialized with one of both.
