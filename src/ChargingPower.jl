module ChargingPower

using SparseArrays

export majorana_operators

include(joinpath(@__DIR__, "building_operators/majorana_operators.jl"))

end # module ChargingPower
