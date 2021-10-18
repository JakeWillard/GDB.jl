
struct AffineMap{T}

    Mat :: SparseMatrixCSC{T, Int64}
    f :: Vector{T}

end


function Base.:*(A::AffineMap{T}, v::Vector{T}) where {T}
    A.Mat*v + A.f
end


function Base.:*(M::SparseMatrixCSC{T, Int64}, A::AffineMap{T}) where {T}
    AffineMap{T}(M*A.Mat, M*A.f)
end

function Base.:*(A::AffineMap{T}, M::SparseMatrixCSC{T, Int64}) where {T}
    AffineMap{T}(A.Mat*M, A.f)
end


function Base.:*(A1::AffineMap{T}, A2::AffineMap{T}) where {T}
    Mat = A1.Mat*A2.Mat
    f = A1.Mat*A2.f + A1.f
    AffineMap{T}(Mat, f)
end


mutable struct LinearProblem{T}

    A :: SparseMatrixCSC{T, Int64}
    Adiag :: Vector{T}
    b :: Vector{T}
    x :: Vector{T}
    r :: Vector{T}

end
function LinearProblem{T}(A, x, b) where {T}
    Adiag = diag(A)
    r = b - A*x
    return LinearProblem{T}(A, Adiag, b, x, r)
end


function compute_residual!(L::LinearProblem{T}) where {T}
    L.r[:] = L.b - L.A*L.x
end


function coarsen_problem!(L::LinearProblem{T}, Restr::AffineMap{T}, Interp::AffineMap{T}) where {T}

    lhs = Restr * (L.A*Interp)
    rhs = Restr * L.b

    A = lhs.Mat
    b = rhs - lhs.f
    x = zeros(length(b))

    return LinearProblem{T}(A, x, b)
end


function jsmooth!(L::LinearProblem{T}, n::Int64) where {T}

    for _=1:n
        L.x[:] += L.r ./ L.Adiag
        compute_residual!(L)
    end
end
