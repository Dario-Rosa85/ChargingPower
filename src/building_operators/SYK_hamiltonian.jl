#This function builds the quadratic SYK Hamiltonian based on a given graph
function SYK_hamiltonian(majorana_operators, graph)
    n_majorana = length(majorana_operators)
    SYK_hamiltonian_matrix = spzeros(2^(n_majorana/2), 2^(n_majorana/2))
    for edge in edges(graph) 
        SYK_hamiltonian_matrix .+= (im/2) * sqrt(n_majorana/(2 * length(edges))) * randn().* majorana_operators[src(edge)] * majorana_operators[dst(edge)]
    end
    return SYK_hamiltonian_matrix
end