#This function builds the quadratic SYK Hamiltonian based on a given graph
function SYK_hamiltonian(majorana_operators, graph)
    n_majorana = length(majorana_operators)
    SYK_hamiltonian_matrix = spzeros(ComplexF64,2^floor(Int, n_majorana/2), 2^floor(Int, n_majorana/2))
    for edge in edges(graph) 
        SYK_hamiltonian_matrix .+= (im/2) * sqrt(n_majorana/(2 * length(edges(graph)))) * randn().* majorana_operators[src(edge)] * majorana_operators[dst(edge)]
        # SYK_hamiltonian_matrix .+= (im/2) * sqrt(1/(n_majorana)) * randn().* majorana_operators[src(edge)] * majorana_operators[dst(edge)]
    end
    return SYK_hamiltonian_matrix
end