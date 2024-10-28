module ChargingPower

using SparseArrays
using Graphs

export majorana_operators, SYK_hamiltonian

include(joinpath(@__DIR__, "building_operators/majorana_operators.jl"))
include(joinpath(@__DIR__, "building_operators/SYK_hamiltonian.jl"))
include(joinpath(@__DIR__, "building_operators/battery_hamiltonians.jl"))

end # module ChargingPower
