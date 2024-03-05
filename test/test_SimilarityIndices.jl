using VegSci
using Test
using NamedArrays

@testset "SimilarityIndices.jl" begin

    x = generate_test_array(rown = 10, coln = 10, meancoloccs = 6, rowprefix = "SiteA-", colprefix = "Species")
    y = generate_test_array(rown = 5, coln = 10, meancoloccs = 6, rowprefix = "SiteB-", colprefix = "Species")
    xy = VegSci.merge_namedarrays([x, y])

    @testset "steinhaus_coefficient" begin
        steinhaus_results = VegSci.steinhaus_coefficient(x, y)
        @test typeof(steinhaus_results) <: NamedMatrix
        @test size(steinhaus_results) == (10, 5)
        @test all(x -> x .>= 0.0, steinhaus_results)
        @test all(x -> x .<= 1.0, steinhaus_results)
        @test names(steinhaus_results)[1] == names(x)[1]
        @test names(steinhaus_results)[2] == names(y)[1]
    end
    @testset "jaccard_coefficient" begin
        jaccard_results = VegSci.jaccard_coefficient(xy)
        @test typeof(jaccard_results) <: NamedMatrix
        @test size(jaccard_results) == (15, 15)
        @test all(x -> x .>= 0.0, jaccard_results)
        @test all(x -> x .<= 1.0, jaccard_results)
        @test names(jaccard_results)[1] == [names(x)[1]; names(y)[1]]
        @test names(jaccard_results)[2] == [names(x)[1]; names(y)[1]]
    end
end