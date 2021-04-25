

abstract type Vorticity <: Physical end
abstract type MagneticPotential <: Physical
abstract type ElectricPotential <: Physical end
abstract type ElectronTemp <: Physical end






function solve(w::Variable{Vorticity}, lin::LinearSystem, bs::Array{Boundary, 1})

end
