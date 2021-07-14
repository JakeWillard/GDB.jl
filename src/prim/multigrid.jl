function transformation_matrices(Nx::Int32, Ny::Int32)
    # calculates the interpolation and restriction matrices given a certain resolution in each direction
    # Nx = grid.Nx
    # Ny = grid.Ny
    coarsexres = ceil((Nx-1)/2)
    coarseyres = ceil((Ny-1)/2)

    xcoarse = zeros(Float64, coarsexres)
    ycoarse = zeros(Float64, coarseyres)  # using the reduction in number of points from the paper
    # way to define these with each new iteration? Another overarching loop with the number of iterations?
    # could do a while loop potentially, if trying to get under a threshold

    # define the interpolation matrix, then the restriction matrix as its transpose
    
    # only working in 1D so far
    Interpx = zeros(Nx, coarsexres)
    Restrx1D = zeros(coarsexres, Nx)
    Restrx = zeros(coarsexres, Nx)

    for i=1:Nx
        for j=1:coarsexres
            # sort of brute force based on the form of the interpolation matrix given in the MIT site
            Interpx[2(j-1)+1, j] = .5
            Interpx[2(j-1)+2, j] = 1
            Interpx[2(j-1)+3, j] = .5
            # may not need outer loop if indexing only depends on j
        end
    end

    Restrx1D = .5*transpose(Interpx) # definition of restriction matrix from the paper
    

    # moving to 2D?
    # defining values based on 2D interpolation scheme in the paper
    """for i=1:Nx
        for j=1:coarsexres
            if isodd(i)
                if isodd(j)
                    Interpx[i,j] = .25
                elseif iseven(j)
                    Interpx[i,j] = .5
                end
            elseif iseven(i)
                if isodd(j)
                    Interpx[i,j] = .5
                elseif iseven(j)
                    Interpx[i,j] = 1
                end
            end
        end
    end"""

    # Restrx = .25*transpose(Interpx) # definition of restriction matrix
    Restr2d = kron(Restrx1D, Restrx1D) # assuming we're using the 1D matrix
    Interp2d = 4*transpose(Restr2D)
    # same Interpolation/Restriction matrices in x and y?

    return Restr2d, Interp2d, coarsexres, coarseyres
end

# not sure if fine_to_coarse and coarse_to_fine are needed, but they could be useful so I'll keep them in
function fine_to_coarse(Nx::Int32, Ny::Int32, x::Vector{Float64}, y::Vector{Float64})
    # restricts values from fine grid to fit them onto coarse grid

    Restriction, Interpolation, newxres, newyres = transformation_matrices(Nx, Ny)

    # calculating the coarse resolution for reducing the number of grid points by 2 in each dimension, assuming odd number of points in each dimension
    coarsexres = ceil((Nx-1)/2)
    coarseyres = ceil((Ny-1)/2)

    # initializing coarse grid values in each dimension
    xcoarse = zeros(Float64, coarsexres)
    ycoarse = zeros(Float64, coarseyres)

    # assuming dimensionality of both are the same? 
    xcoarse[:] = Restriction * x
    ycoarse[:] = Restriction * y

    return xcoarse, ycoarse
end

function fine_to_coarse(Nx::Int32, Ny::Int32, x::Vector{Float64}, y::Vector{Float64})
    # projects vector values from the fine grid to the coarse grid using interpolation

    Restriction, Interpolation, newxres, newyres = transformation_matrices(Nx, Ny)

    # calculating the coarse resolution for reducing the number of grid points by 2 in each dimension, assuming odd number of points in each dimension
    coarsexres = ceil((Nx-1)/2)
    coarseyres = ceil((Ny-1)/2)

    # initializing coarse grid values in each dimension
    xcoarse = zeros(Float64, coarsexres)
    ycoarse = zeros(Float64, coarseyres)

    # assuming dimensionality of both are the same? 
    xcoarse[:] = Interpolation * x
    ycoarse[:] = Interpolation * y

    return xcoarse, ycoarse
end

function project_to_coarse(A::SparseMatrixCSC, Nx::Int32, Ny::Int32, R::Matrix{Float64}, I::Matrix{Float64})
    # projects the matrix from the matrix equation Ax = b onto the coarse grid

    # acquiring restriction and interpolation matrices
    # R, I, xres, yres = transformation_matrices(Nx, Ny)

    # using formula for projection of matrix onto coarse grid
    Acoarse = R * A * I

    return Acoarse
end

function correction_matrix(A::SparseMatrixCSC}, Nx::Int32, Ny::Int32)
    # calculates the error correction matrix s

    # acquiring restriction and interpolation matrices
    R, I, xres, yres = transformation_matrices(Nx, Ny)

    # carrying out explicit calculation of s
    S = I * inv(R * A * I) * R * A

    return S
end

# NOTE THAT THE x I AM USING HERE IS THE ACTUAL SOLUTION AND NOT THE GUESSED/KNOWN SOLUTION
# Also the x above is x as in the coordinate, below the x is the solution to the matrix equation Ax = b

function residual(A::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64})
    # calculates the residual associated with the vector equation we're trying to solve

    # definition of residual on fine grid
    rfine = b - A * x

    return rfine
end

function restricted_residual(A::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64}, Nx::Int32, Ny::Int32)
    # calculates the residual on the coarse grid

    # getting restriction/interpolation matrices on the given grid
    R, I, xres, yres = transformation_matrices(Nx::Int32, Ny::Int32)

    # obtaining residual on fine grid
    rfine = residual(A, x, b)

    # calculation of restriction residual
    rcoarse = R * rfine

    return rcoarse
end

function multigrid_error(e::Vector{Float64}, A::SparseMatrixCSC}, Nx::Int32, Ny::Int32)
    # given the error e = x(guessed) - x(calculated), this gives the multigrid error correction

    # calculation of correction matrix S that maps from coarse to fine grid
    S = correction_matrix(A, Nx, Ny)

    # definition of multigrid correction given true fine grid error
    E = S * e
    
    return E
end

# Need to solve Acoarse = Ecoarse * rcoarse, then interpolate back and add Efine to x(calculated)

# may not need y vector - check later if issue doesn't get resolved
function multigrid_step(grid::Grid, A::SparseMatrixCSC, xcalc::Vector{Float64}, xguess::Vector{Float64}, b::Vector{Float64})
    # steps through one step of multigrid from fine to coarse and back according to the steps in the MIT page
    
    Nx = grid.Nx
    Ny = grid.Ny

    R1, I1, Nx1, Ny1 = transformation_matrices(Nx, Ny)
    R2, I2, Nx2, Ny2 = transformation_matrices(Nx1, Ny1)
    R3, I3, Nx3, Ny3 = transformation_matrices(Nx2, Ny2)

    rcoarse1 = restricted_residual(A, xcalc, b, Nx, Ny)
    Acoarse1 = project_to_coarse(A, Nx, Ny, R1, I1)
    xcalc1 = R1 * xcalc
    b1 = R1 * b

    rcoarse2 = restricted_residual(Acoarse1, xcalc1, b1, Nx1, Ny1)
    Acoarse2 = project_to_coarse(Acoarse1, Nx1, Ny1, R2, I2)
    xcalc2 = R2 * xcalc
    b2 = R2 * b

    rcoarse3 = restricted_residual(Acoarse2, xcalc2, b2, Nx2, Ny2)
    Acoarse3 = project_to_coarse(Acoarse2, Nx2, Ny2, R3, I3)
    xcalc3 = R3 * xcalc
    b3 = R3 * b


    # actual calculation - use a vcycle for example?
    # iterate at each step, carry over E into next grid and use it as guess for iteration on that grid
    # potentially start with exact solution at lowest grid level
    # for example
    Ecoarse3 = Acoarse3 \ rcoarse3
    Ecoarse2 = I3 * Ecoarse3
    Ecoarse2_1 = p_jacobi(Acoarse2, Ecoarse2, rcoarse2, 1, 1, 1, .1)
    Ecoarse1 = I2 * Ecoarse2_1
    Ecoarse1_1 = p_jacobi(Acoarse1, Ecoarse1, rcoarse1, 1, 1, 1, .1)
    Efine = I1 * Ecoarse1_1

    xcalc = xcalc + Efine

    return xcalc # may want to return more depending on the information we want
end


function v_cycle(grid::Grid, A::SparseMatrixCSC, xcalc::Vector{Float64}, b::Vector{Float64})
    # steps through one step of multigrid from fine to coarse and back according to the steps in the MIT page
    
    Nx = grid.Nx
    Ny = grid.Ny

    R1, I1, Nx1, Ny1 = transformation_matrices(Nx, Ny)
    R2, I2, Nx2, Ny2 = transformation_matrices(Nx1, Ny1)
    R3, I3, Nx3, Ny3 = transformation_matrices(Nx2, Ny2)

    rcoarse1 = restricted_residual(A, xcalc, b, Nx, Ny)
    Acoarse1 = project_to_coarse(A, Nx, Ny, R1, I1)
    xcalc1 = R1 * xcalc
    b1 = R1 * b

    rcoarse2 = restricted_residual(Acoarse1, xcalc1, b1, Nx1, Ny1)
    Acoarse2 = project_to_coarse(Acoarse1, Nx1, Ny1, R2, I2)
    xcalc2 = R2 * xcalc1
    b2 = R2 * b1

    rcoarse3 = restricted_residual(Acoarse2, xcalc2, b2, Nx2, Ny2)
    Acoarse3 = project_to_coarse(Acoarse2, Nx2, Ny2, R3, I3)
    xcalc3 = R3 * xcalc2
    b3 = R3 * b2


    Einitial = zeros({Float64}, size(xcalc))
    Ecoarse1 = jacobi_smooth(Acoarse1, Einitial, rcoarse1, 1, 1, 1, 1, .1)
    Ecoarse2 = R2 * Ecoarse1
    Ecoarse2_1 = jacobi_smooth(Acoarse2, Ecoarse2, rcoarse2, 1, 1, 1, 1, .1)
    Ecoarse3 = R3 * Ecoarse2_1
    Ecoarse3_1 = jacobi_smooth(Acoarse3, Ecoarse3, rcoarse3, 1, 1, 1, 1, .1)
    Ecoarse2_2 = I3 * Ecoarse3_1
    Ecoarse2_3 = jacobi_smooth(Acoarse2, Ecoarse2_2, rcoarse2, 1, 1, 1, 1, .1)
    Ecoarse1_2 = I2 * Ecoarse2_3
    Ecoarse1_3 = jacobi_smooth(Acoarse1, Ecoarse1_2, rcoarse1, 1, 1, 1, 1, .1)
    Efine = I1 * Ecoarse1_3


    xcalc = xcalc + Efine

    return xcalc # may want to return more depending on the information we want
end