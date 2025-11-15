using VegSci
using Test
using NamedArrays

@testset "Clustering.jl" begin
    ta = VegSci.generate_test_array(rown = 15, coln = 30, meancoloccs = 10, rowprefix = "SiteA-", colprefix = "Species")
    nm = VegSci.prop_to_domin(ta)'
    @testset "hamming_distance" begin
        hd_results = VegSci.hamming_distance(nm)
        @test typeof(hd_results) <: NamedMatrix
        @test size(hd_results) == (15, 15)
        @test all(x -> x >= 0, hd_results)
        @test names(hd_results)[1] == names(nm)[2]
        @test names(hd_results)[2] == names(nm)[2]
    end
end