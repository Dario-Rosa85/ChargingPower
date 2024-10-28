using ChargingPower, Test, LinearAlgebra

@testset "majorana_algebra" begin
    n_matrices = 8
    matrices = majorana_operators(n_matrices)
    for i in 1:n_matrices  
        for j in 1:n_matrices 
            if i != j
                @test isapprox(matrices[i] * matrices[j] + matrices[j] * matrices[i], 0.0 .* Matrix(I, Int64(2^(n_matrices/2)), Int(2^(n_matrices/2))))
            else
                @test isapprox(matrices[i] * matrices[j] + matrices[j] * matrices[i], 2.0 .* Matrix(I, Int64(2^(n_matrices/2)), Int(2^(n_matrices/2))))
            end 
        end
    end
end