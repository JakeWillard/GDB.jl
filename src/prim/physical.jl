
abstract type Physical end

mutable struct Variable{T<:Physical}

    value::Array{Float64, 3}
    x::Matrix{Float64}
    xx::Matrix{Float64}
    xy::Matrix{Float64}
    y::Matrix{Float64}
    yy::Matrix{Float64}
    t::Matrix{Float64}

end


function Variable{T}(v0::Vector{Float64}, Nz::Int) where {T<:Physical}

    Nk = length(v0)
    value = zeros(Float64, (Nk, Nz, 3))
    for i=1:Nz
        for j=1:3
            value[:,i,j] = v0[:]
        end
    end

    blank = zeros(Float64, (Nk, Nz))
    return Variable{T}(value, blank, blank, blank, blank, blank, blank)
end


function Variable{T}(f::Function, Nz::Int, grid::Grid) where {T<:Physical}

    Nk = grid.Nk
    v0 = zeros(Float64, Nk)

    for k=1:Nk
        x = grid.points[1,k]
        y = grid.points[2,k]
        v0[k] = f(x, y)
    end

    return Variable{T}(v0, Nz)
end


function to_readable(var::Variable, grid::Grid)

    vals = grid.inverse_projection * var.value[:,:,2]

    Nx = grid.Nx
    Ny = grid.Ny
    Nz = size(vals)[2]
    out = zeros(Float64, (grid.Nx, grid.Ny, Nz))

    for j=1:Ny
        out[:,j,:] = vals[1+(j-1)*Nx:j*Nx,:]
    end

    return out
end
