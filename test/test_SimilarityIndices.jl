using VegSci
using Test
using NamedArrays

@testset "SimilarityIndices.jl" begin

    x = generate_test_array(rown = 10, coln = 10, meancoloccs = 6, rowprefix = "SiteA-", colprefix = "Species")
    y = generate_test_array(rown = 5, coln = 10, meancoloccs = 6, rowprefix = "SiteB-", colprefix = "Species")

    @testset "czekanowski_index" begin
        czeksim_results = VegSci.czekanowski_index(x, y)
        @test typeof(czeksim_results) <: NamedMatrix
        @test size(czeksim_results) == (10, 5)
        @test all(x -> x .>= 0.0, czeksim_results)
        @test all(x -> x .<= 1.0, czeksim_results)
        @test names(czeksim_results)[1] == names(x)[1]
        @test names(czeksim_results)[2] == names(y)[1]
    end
end