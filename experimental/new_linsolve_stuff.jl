
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
    Base.:*(l::LinearMap, v::Vector{Float64}) = l.M * v + l.f
    LinearMap(A::SparseMatrixCSC, b::Vector{Float64}) = new(A, b)
end


function test_jacobi(A::SparseMatrixCSC, x0, b, w, Nout, Nin, nchunks)

    if nchunks > length(workers())
        nchunks = 1
    end

    D = w*Diagonal(1 ./ diag(A))
    J = LinearMap(I - D*A, Vector(D*b))
    x = x0[:]

    for _=1:Nout

        xnext = DArray((length(x0),), workers()[1:nchunks], nchunks) do inds
            Jp = J[inds[1]]
            xp = x[:]
            for __=1:Nin
                xp[inds[1]] = Jp * xp
            end
            xp[inds[1]]
        end
        x[:] = Array(xnext)
        @info err = norm(A*x - b)
    end

    return x
end




function require_boundary_conditions(L::LinearMap, xb::Vector{Float64}, gc::GhostConditions)

    Anew, Bterm = require_boundary_conditions(L.M, gc)
    fnew = gc.Proj*L.f + Bterm*xb

    return LinearMap(Anew, fnew)
end


function solve_pde(A::SparseMatrixCSC, b::Vector{Float64}, xb::Vector{Float64}, gc::GhostConditions)

    L = require_boundary_conditions(LinearMap(A, -b), xb, gc)
    x = L.M \ -L.f

    return extrapolate_ghosts(x, xb, gc)
end
