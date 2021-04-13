

abstract type ModelVariable end
abstract type Temp <: ModelVariable end
abstract type Density <: ModelVariable end

struct Variable{T<:ModelVariable}

    value::Matrix{Float64}
    x::Matrix{Float64}
    xx::Matrix{Float64}
    xy::Matrix{Float64}
    y::Matrix{Float64}
    yy::Matrix{Float64}
    t::Matrix{Float64}

end


function Variable{T}(value::Matrix{Float64}) where {T}
    blank = zeros(Float64, size(value))
    return Variable{T}(value, blank, blank, blank, blank, blank, blank)
end
