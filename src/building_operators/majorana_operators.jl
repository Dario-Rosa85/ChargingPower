# This function build Majorana operators satisfying the relation \left\{\gamma_i, \gamma_j\right\} = 2 \delta_{ij}
function majorana_operators(n_majorana)
    pauli_x = sparse([0 1; 1 0])
    pauli_y = sparse([0 -1 * im; 1 * im 0])
    pauli_z = sparse([1 0; 0 -1])
    identity = sparse([1 0; 0 1])
    majorana_matrices = SparseMatrixCSC{ComplexF64,Int64}[]
    for i in 1:n_majorana 
        if i == 1
            push!(majorana_matrices, pauli_x)
        elseif i == 2
            push!(majorana_matrices, pauli_y)
        else
            push!(majorana_matrices, pauli_z)
        end
    end
    @inbounds Threads.@threads :static for i in 1:2:n_majorana-1 
        for j in 2:floor(Int, n_majorana/2)
            if 2*j - 1 < i
                majorana_matrices[i] = kron(majorana_matrices[i], pauli_z)
            elseif 2*j - 1 == i
                majorana_matrices[i] = kron(majorana_matrices[i], pauli_x)
            else
                majorana_matrices[i] = kron(majorana_matrices[i], identity)
            end             
        end
    end
    @inbounds Threads.@threads :static for i in 2:2:n_majorana 
        for j in 2:floor(Int, n_majorana/2)
            if 2*j < i
                majorana_matrices[i] = kron(majorana_matrices[i], pauli_z)
            elseif 2*j == i
                majorana_matrices[i] = kron(majorana_matrices[i], pauli_y)
            else
                majorana_matrices[i] = kron(majorana_matrices[i], identity)
            end             
        end
    end
    return majorana_matrices    
end