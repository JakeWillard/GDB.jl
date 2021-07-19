
mutable struct MultiResMatrix <: AbstractMatrix{Float64}

    mats :: Array{SparseMatrixCSC{Float64, Int32}, 1}
    lvl :: Int32

    MultiResMatrix(mats) = MultiResMatrix(mats, 1)
    Base.size(m::MultiResMatrix) = Base.size(m.mats[m.level])
    Base.getindex(m::MultiResMatrix, inds...) = Base.getindex(m.mats[m.level], inds...)
    Base.setindex!(m::MultiResMatrix, X, inds...) = Base.setindex!(m.mats[m.level], X, inds...)
    uplvl!(m::MultiResMatrix) = begin m.lvl += 1 end
    dwnlvl!(m::MultiResMatrix) = begin m.lvl -= 1 end
end


function MultiResMatrix(A::SparseMatrixCSC, grd::Grid)

    As = [A]
    for l=1:length(grd.restrictions)
        Restr = grd.restrictions[l]
        Interp = grd.interpolations[end-(l-1)]
        append!(As, Restr * As[end] * Interp)
    end

    return MultiResMatrix(As)
end


function MultiResMatrix(f::Function, A::SparseMatrixCSC, grd::Grid)

    Fs = [f(A)]
    for l=1:length(grd.restrictions)
        Restr = grd.restrictions[l]
        Interp = grd.interpolations[end-(l-1)]
        append!(Fs, f(Restr * As[end] * Interp))
    end

    return MultiResMatrix(Fs)
end


function multires_jacobi_smoother(A::SparseMatrixCSC, grd::Grid, w::Float64)

    f(X) = I - w*inv(Diagonal(X)) * X
    g(X) = inv(Diagonal(X))

    M = MultiResMatrix(f, A, grd)
    Dinv = MultiResMatrix(g, A, grd)
    return M, Dinv
end


function jacobi_vcycle(M::MultiResMatrix, Dinv::MultiResMatrix, )
