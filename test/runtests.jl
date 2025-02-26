using RandomizedSparsification
using Test
using Aqua

@testset "RandomizedSparsification.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(RandomizedSparsification; piracies=(; broken=true))
    end
    @testset "rank" begin
        using LinearAlgebra
        using RandomizedSparsification: srank, rankrand
        A_large = rankrand(1000, 1500, 20)
        @test rank(A_large) == 20
        @test srank(A_large) <= rank(A_large)
    end
    @testset "FindIndices" begin
        using RandomizedSparsification: FindIndices as FI

        P_ij = [1 2 3; 4 5 6; 7 8 9] / sum(1:9)
        P_ij_cs = cumsum(view(P_ij, :))
        rand01_vec = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        out = [
            CartesianIndex(1, 1), # 0.01 <= 0.0222222
            CartesianIndex(2, 1), # 0.0222222 < 0.1 <= 0.1111111 
            CartesianIndex(3, 1), # 0.1111111 < 0.2 <= 0.266667  
            CartesianIndex(1, 2), # 0.266667 < 0.3 <= 0.311111
            CartesianIndex(2, 2), # 0.311111 < 0.4 <= 0.422222 
            CartesianIndex(3, 2), # 0.422222 < 0.5 <= 0.6
            CartesianIndex(3, 2), # 0.422222 < 0.6 <= 0.6
            CartesianIndex(2, 3), # 0.666667 < 0.7 <= 0.8
            CartesianIndex(2, 3), # 0.666667 < 0.8 <= 0.8
            CartesianIndex(3, 3), # 0.8 < 0.9
        ]

        @test FI._indices_sequencial(P_ij_cs, rand01_vec, 3, 3) == out
        @test FI._indices_groupsort(P_ij_cs, rand01_vec, 3, 3) == out
    end
end
