module ChargingPower

using SparseArrays
using Graphs
using KrylovKit

export majorana_operators, SYK_hamiltonian, Z_model_hamiltonian, X_model_hamiltonian, fake_X_model_hamiltonian

include(joinpath(@__DIR__, "building_operators/majorana_operators.jl"))
include(joinpath(@__DIR__, "building_operators/SYK_hamiltonian.jl"))
include(joinpath(@__DIR__, "building_operators/battery_hamiltonians.jl"))

end # module ChargingPower
