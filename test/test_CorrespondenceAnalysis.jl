using VegSci
using Test
using NamedArrays

@testset "CorrespondenceAnalysis.jl" begin

    x = VegSci.generate_test_array(rown = 30, coln = 10, meancoloccs = 7, rowprefix = "", colprefix = "Species")

    @testset "correspondence_analysis" begin
        ca_results = VegSci.correspondence_analysis(x)
        @test typeof(ca_results) <: NamedTuple
        @test keys(ca_results) == (:sv, :rownames, :rowmass, :rowcoord, :colnames, :colmass, :colcoord, :N)
    end
end