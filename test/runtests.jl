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
end
