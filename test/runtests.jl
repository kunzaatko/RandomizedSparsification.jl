using RandomizedSparsification
using Test
using Aqua

@testset "RandomizedSparsification.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(RandomizedSparsification)
    end
    # Write your tests here.
end
