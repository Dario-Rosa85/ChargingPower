#This function builds the quadratic SYK Hamiltonian based on a given graph
function SYK_hamiltonian(majorana_operators, graph)
    n_majorana = length(majorana_operators)
    SYK_hamiltonian_matrix = spzeros(ComplexF64,2^floor(Int, n_majorana/2), 2^floor(Int, n_majorana/2))
    if length(edges(graph)) == n_majorana * (n_majorana - 1) /2
        normalization_factor = (im/2) * sqrt(1/(n_majorana))
    else
        normalization_factor = (im/2) * sqrt(n_majorana/(2 * length(edges(graph))))
    end
    for edge in edges(graph) 
        SYK_hamiltonian_matrix .+= normalization_factor * randn().* majorana_operators[src(edge)] * majorana_operators[dst(edge)]
    end
    return SYK_hamiltonian_matrix
end

#This function produces the SYK Hamiltonian on a twisted circle
function twisted_circle_SYK_hamiltonian(majorana_operators, graph)
    n_majorana = length(majorana_operators)
    SYK_hamiltonian_matrix = spzeros(ComplexF64,2^floor(Int, n_majorana/2), 2^floor(Int, n_majorana/2))
    normalization_factor = (im/2) * sqrt(n_majorana/(2 * length(edges(graph))))
    for edge in edges(graph)
        if src(edge) + dst(edge) == n_majorana + 1
            SYK_hamiltonian_matrix .+= normalization_factor * randn().* majorana_operators[src(edge)] * majorana_operators[dst(edge)]
        elseif src(edge) + dst(edge) < n_majorana + 1
            to_be_changed = max(src(edge), dst(edge))
            to_remain = min(src(edge), dst(edge))
            SYK_hamiltonian_matrix .+= normalization_factor * randn().* majorana_operators[to_remain] * majorana_operators[n_majorana + 1 - to_be_changed]
        else
            to_be_changed = min(src(edge), dst(edge))
            to_remain = max(src(edge), dst(edge))
            SYK_hamiltonian_matrix .+= normalization_factor * randn().* majorana_operators[n_majorana + 1 - to_be_changed] * majorana_operators[to_remain]
        end 
    end
    return SYK_hamiltonian_matrix
end