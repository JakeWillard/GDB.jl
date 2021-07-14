
macro variable2d(name)

    type_def = Meta.parse("""
    mutable struct $(name){T, N} <: AbstractArray{T, N}
        vals :: Array{T, N}

        Base.size(v::$(name)) = Base.size(v.vals)
        Base.getindex(v::$(name), inds...) = Base.getindex(v.vals, inds...)
        Base.setindex!(v::$(name), X, inds...) = Base.setindex!(v.vals, X, inds...)
        $(name)(x::Vector{Float64}) = $(name){Float64, 2}(hcat(x, x))
    end
    """)

    return type_def
end
