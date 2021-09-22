### ImmersedMirrors.jl Guide

## Dependencies

* ProgressMeter.jl
* ForwardDiff.jl
* Logging.jl
* RecipesBase.jl

## Overview

The purpose of this module is to make discretizing problems on a square lattice with arbitrary immersed boundary as simple as it would be for a rectangular boundary. Any linear operators that you can define on a rectangular grid can be easily "translated" to a grid with trimmed edges, which cuts down on the extra computational load associated with immersed boundary problems that keep a rectangular grid. Linear operators that correspond to mirroring values across an immersed boundary can be easily calculated, which are useful in imposing symmetry and anti-symmetry boundary conditions.

## Grid Constructor

A Grid struct is initialized given three pieces of info: (1) bounding box for the gridpoint search, (2) grid resolution, and (3) rule for determining what points within the bounding box exist within the domain. You specify this info with the following inputs:
* inside -- a function of (x,y) which returns true for points inside the domain and false for points outside the domain.
* r0 -- vector for the bottom left corner of the bounding box.
* r1 -- vector for the top right corner of the bounding box.
* dr -- distance between neighboring gridpoints (the lattice is square, not just rectangular, in this newest version).
* Nz -- number of planes in the phi direction.

The Grid constructor first uses r0, r1, and dr to determine the number of gridpoints in each direction Nx and Ny to fill the bounding box. This means that the input r1 is actually treated as a guess/approximation, and the true r1 is actually determined to be the number closest to input r1 such that Ny = r1[2]/dr is an integer. Next, the constructor loops through each point in the bounding box and checks which ones are inside the domain. Which of the points exist within the domain, and where they show up in the list, is used to set each of the attributes of the Grid struct.

To initialize a Grid in the code, you can simply call Grid(inside, r0, r1, dr, Nz), or use the do-block syntax:

    grd = Grid(dr, Nz) do
      # code that defines r0, r1, and inside
      inside, r0, r1
    end

You can also of course just use the do-block to define the function inside:

    grd = Grid(r0, r1, dr, Nz) do x,y
       # function that returns true when (x,y) is inside the domain
    end

## Functions on the Grid

## Operators on the Grid

## Mirror Constructor

## GhostData Constructor

## Setting values at Ghost points.

## Imposing Boundary Conditions

## Examples
