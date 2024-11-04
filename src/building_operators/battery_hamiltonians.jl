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
            hamiltonian_matrix .+=  majorana_operators[1]
        else
            provisional_matrix =  majorana_operators[1]
            for j in 2:(2 * k - 1) 
                provisional_matrix *=  majorana_operators[j]
            end
            hamiltonian_matrix .+=  (-im)^(k + 1) .* provisional_matrix
        end
    end
    return hamiltonian_matrix
end

# routine used to renormalize an Hamiltonian to lie between -1 and 1
function renormalization_hamiltonian!(hamiltonian_to_renormalize)
    matrix_size = size(hamiltonian_to_renormalize, 2)
    sparse_identity = spzeros(Float64, matrix_size, matrix_size)
    @inbounds for i in 1:matrix_size
        sparse_identity[i, i] = 1.0
    end
    e_max = real(first(first(eigsolve(hamiltonian_to_renormalize, rand(ComplexF64, matrix_size), 1, :LR))))
    e_min = real(first(first(eigsolve(hamiltonian_to_renormalize, rand(ComplexF64, matrix_size), 1, :SR))))
    rescaling_factor = 2 / (e_max - e_min)
    shifting_factor = (e_max + e_min) / (e_max - e_min)
    hamiltonian_to_renormalize .= rescaling_factor .* hamiltonian_to_renormalize .- shifting_factor .* sparse_identity
end

# This function creates the Hamiltonian for the fake X-model
function fake_X_model_hamiltonian(majorana_operators)
    n_majorana = length(majorana_operators)
    hamiltonian_matrix = spzeros(ComplexF64, 2^floor(Int, n_majorana/2), 2^floor(Int, n_majorana/2))
    for k in 1: floor(Int, n_majorana/2) 
        if k == 1
            hamiltonian_matrix .+= majorana_operators[1]
        else
            provisional_matrix =  majorana_operators[1]
            for j in 2:k 
                provisional_matrix *=  majorana_operators[2*j - 1]
            end
            hamiltonian_matrix .+= (-im)^(k + 1) .* provisional_matrix
        end
    end
    return n_majorana/2 .* renormalization_hamiltonian!(hamiltonian_matrix)
end