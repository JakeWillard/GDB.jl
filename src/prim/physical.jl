
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

    Nxy = length(v0)
    value = zeros(Float64, (Nxy, Nz, 3))
    for i=1:Nz
        for j=1:3
            value[:,i,j] = v0[:]
        end
    end

    blank = zeros(Float64, (Nxy, Nz))
    return Variable{T}(value, blank, blank, blank, blank, blank, blank)
end



function Variable{T}(f::Function, grid::Grid, Nz::Int) where {T<:Physical}

    Nxy = grid.Nk
    v0 = zeros(Float64, grid.Nk)
    for k=1:grid.Nk
        x = grid.points[1,k]
        y = grid.points[2,k]
        v0[k] = f(x, y)
    end

    return Variable{T}(v0, Nz)
end
