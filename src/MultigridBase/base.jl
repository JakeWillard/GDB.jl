
struct AffineMap

    Mat :: SparseMatrixCSC{Float64, Int64}
    f :: Vector{Float64}

end


function Base.:*(A::AffineMap, v::Vector)
    A.Mat*v + A.f
end


function Base.:*(M::SparseMatrixCSC{Float64, Int64}, A::AffineMap)
    AffineMap(M*A.Mat, M*A.f)
end

function Base.:*(A::AffineMap, M::SparseMatrixCSC{Float64, Int64})
    AffineMap(A.Mat*M, A.f)
end


function Base.:*(A1::AffineMap, A2::AffineMap)
    Mat = A1.Mat*A2.Mat
    f = A1.Mat*A2.f + A1.f
    AffineMap(Mat, f)
end


mutable struct LinearProblem

    A :: SparseMatrixCSC{Float64, Int64}
    Adiag :: Vector{Float64}
    b :: Vector{Float64}
    x :: Vector{Float64}
    r :: Vector{Float64}

end
function LinearProblem(A::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64})
    Adiag = diag(A)
    r = b - A*x
    return LinearProblem(A, Adiag, b, x, r)
end


function compute_residual!(L::LinearProblem)
    L.r[:] = L.b - L.A*L.x
end


function coarsen_problem!(L::LinearProblem, Restr::AffineMap, Interp::AffineMap)

    lhs = Restr * (L.A*Interp)
    rhs = Restr * L.b

    A = lhs.Mat
    b = rhs - lhs.f
    x = zeros(length(b))

    return LinearProblem(A, x, b)
end


function jsmooth!(L::LinearProblem, n::Int64)

    for _=1:n
        L.x[:] += L.r ./ L.Adiag
        compute_residual!(L)
    end
end
