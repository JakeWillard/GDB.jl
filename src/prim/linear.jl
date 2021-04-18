

struct LinearSystem{T<:Physical}

    A::SparseMatrixCSC
    b::Vector{Float64}

end


function swap_b(lin::LinearSystem, b::Vector{Float64})

    A = lin.A
    return LinearSystem(A, b)
end
