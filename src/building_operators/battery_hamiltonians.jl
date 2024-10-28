# This function creates the Hamiltonian for the Z-model
function Z_model_hamiltonian(majorana_operators)
    n_majorana = length(majorana_operators)
    hamiltonian_matrix = spzeros(ComplexF64,2^floor(Int, n_majorana/2), 2^floor(Int, n_majorana/2))
    for i in 1:floor(Int, n_majorana/2) 
        hamiltonian_matrix .+= -1 * im .* majorana_operators[2 * i - 1] * majorana_operators[2 * i] 
    end
    return hamiltonian_matrix
end

# This function creates the Hamiltonian for the X-model
function X_model_hamiltonian(majorana_operators)
    n_majorana = length(majorana_operators)
    hamiltonian_matrix = spzeros(ComplexF64, 2^floor(Int, n_majorana/2), 2^floor(Int, n_majorana/2))
    for k in 1: floor(Int, n_majorana/2) 
        if k == 1
            hamiltonian_matrix .+= (1/sqrt(2)) .* majorana_operators[1]
        else
            provisional_matrix = sqrt(2) .* majorana_operators[1]
            for j in 2:(2 * k - 1) 
                provisional_matrix *= ((1/sqrt(2)) .* majorana_operators[j])
            end
            hamiltonian_matrix .+= 2^(k - 1/2) * (-1)^(k + 1) .* provisional_matrix
        end
    end
    return hamiltonian_matrix
end