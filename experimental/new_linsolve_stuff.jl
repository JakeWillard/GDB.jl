

struct JacobiSmoother

    C :: SparseMatrixCSC
    f :: Vector{Float64}
    asynch_steps :: Int64
    slices :: Vector{UnitRange}

end


function JacobiSmoother(A::SparseMatrixCSC, b::Vector{Float64}, w::Float64, pchunks, Na, grd::Grid)

    D = w*Diagonal(1 ./ diag(A))
    C = I - D*A
    f = D * b

    pslices = collect(Interators.partition(1:grd.Nk, pchunks))
    slices = UnitRange[]
    for z=1:grd.Nz
        slices = [slices; pslices .+ (z-1)*grd.Nk]
    end

    return JacobiSmoother(C, f, Na, slices)
end


function Base.:*(J::JacobiSmoother, v::Vector{Float64})

    out = pmap(J.slices) do inds
        x = v[:]
        C_loc = J.C[inds, :]
        f_loc = J.f[inds]

        for _=1:J.asynch_steps
            x[inds] = C_loc * x + f_loc
        end
        x[inds]
    end

    return vcat(out)
end


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
