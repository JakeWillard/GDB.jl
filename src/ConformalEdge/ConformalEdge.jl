

module ConformalEdge

export DMap, trace_fieldline, SOL, fluxmap

using LinearAlgebra, ForwardDiff, ConformalMaps, FFTW

include("./dmap.jl")
include("./fieldlines.jl")
include("./sol.jl")

end
