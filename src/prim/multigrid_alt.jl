
# the idea behind this struct is that instances represent a matrix at multiple resolutions. Not only
# can the same matrix at different resolutions be represented by the same variable, but the algebra can be
# simplified in a more elegant way too.
struct MultiResMatrix

    mats :: Array{SparseMatrixCSC{Float64, Int32}, 1}
    lvl_by_row_size :: Dict{Int64, Int64}
    lvl_by_column_size :: Dict{Int64, Int64}

end

# no need for Julia to see this struct as a custom array, but we can overload the multiplication
# function so that multiplying from the left or the right automatically works. Both M*X and X*M will multiply
# X by the correct version of M.
function Base.:*(A::MultiResMatrix, X)
    lvl = A.lvl_by_row_size[size(X)[1]]
    return A.mats[lvl] * X
end
function Base.:*(X, A::MultiResMatrix)
    lvl = A.lvl_by_column_size[size(X)[2]]
    return X * A.mats[lvl]
end


function MultiResMatrix(A::SparseMatrixCSC, grd::Grid)

    As = SparseMatrixCSC[A]
    lvl_by_row_size = Dict(size(A)[2] => 1)
    lvl_by_column_size = Dict(size(A)[1] => 1)


    for l=1:length(grd.restrictions)
        Restr = grd.restrictions[l]
        Interp = grd.interpolations[end-(l-1)]
        Ap = Restr * As[end] * Interp
        append!(As, [Ap])

        lvl_by_row_size[size(Ap)[2]] = l+1
        lvl_by_column_size[size(Ap)[1]] = l+1
    end

    return MultiResMatrix(As, lvl_by_row_size, lvl_by_column_size)
end


function setup_jacobi_smoothing(A::SparseMatrixCSC, w::Float64, grd::Grid)

    A_mr = MultiResMatrix(A, grd)
    Ms = SparseMatrixCSC[]
    Ds = SparseMatrixCSC[]

    for a in A_mr.mats
        D = w*inv(Diagonal(a))
        M = I - w*D*a
        append!(Ms, [M])
        append!(Ds, [D])
    end

    M_mr = MultiResMatrix(Ms, A_mr.lvl_by_row_size, A_mr.lvl_by_column_size)
    D_mr = MultiResMatrix(Ds, A_mr.lvl_by_row_size, A_mr.lvl_by_column_size)

    return A_mr, M_mr, D_mr
end


function jacobi_smooth(n::Int64, M::MultiResMatrix, D::MultiResMatrix, r::Vector{Float64})

    x = zeros(Float64, length(r))
    f = D * r
    for dummy=1:n
        x = M * x + f
    end

    return x
end



function residual(A::MultiResMatrix, x::Vector{Float64}, b::Vector{Float64})
    return b - A * x
end



function vcycle(A::MultiResMatrix, x0::Vector{Float64}, b::Vector{Float64}, M::MultiResMatrix, D::MultiResMatrix, n::Int64, grd::Grid)

    # compute initial residual
    r = residual(A, x0, b)

    # relax on coarser grids
    for Restr in grd.restrictions
        x = jacobi_smooth(n, M, D, r)
        r = Restr*residual(A, x, r)
    end

    # relax on finer grids
    for Interp in grd.interpolations
        x = jacobi_smooth(n, M, D, r)
        r = Interp*residual(A, x, r)
    end

    # final smoothing
    x = jacobi_smooth(n, M, D, r)
    return x0 + x
end


function multigrid_solve(A0::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64}, w::Float64, grd::Grid)

    # setup multiresolution matrices
    A, M, D = setup_jacobi_smoothing(A0, w, grd)

    # compute norm of b, initialize error
    bnorm = norm(b)
    err = 1.0
    println()

    # interate until norm(residual) / bnorm < 10^-8
    while err > 1e-8
        x = vcycle(A, x, b, M, D, 100, grd)
        err = norm(residual(A, x, b)) / bnorm
        println(err)
    end

    return x
end
