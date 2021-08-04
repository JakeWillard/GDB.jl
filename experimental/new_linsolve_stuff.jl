
"""
changes:
    - reflections should be calculated exclusively on ghost points.
    - reflection_matrix() function should have input argument determining if we want
      a reflection or anti-reflection operation, and also input for what to do on other points. 
"""


struct LinearMap

    M :: SparseMatrixCSC
    f :: Vector{Float64}

    Base.size(l::LinearMap) = length(l.f)
    Base.getindex(l::LinearMap, inds) = LinearMap(l.M[inds,:], l.f[inds])
    Base.:*(l::LinearMap, v::Vector{Float64}) = l.M * v + f
    LinearMap(A::SparseMatrixCSC, b::Vector{Float64}) = new(A, b)
end
opnorm(l::LinearMap) = opnorm(Matrix(l.M))   # WARNING: never call this unless the matrix is of small to moderate size.



function subspace_projection(isnt_ghost::Function, grd::Grid)

    is = Int64[1:grd.Nk...]
    js = zeros(Int64, grd.Nk)
    dat = ones(Int64, grd.Nk)
    new_points = zeros(2, grd.Nk)

    k = 0
    for i=1:grd.Nk
        if isnt_ghost(grd.points[:,i]...)
            k += 1
            js[k] = i
            new_points[:,k] = grd.points[:,i]
        end
    end

    return new_points[:,1:k], sparse(is[1:k], js[1:k], dat[1:k], k, grd.Nk)
end
