using VegSci
using Test
using NamedArrays

@testset "CorrespondenceAnalysis.jl" begin

    x = VegSci.generate_test_array(rown = 30, coln = 10, meancoloccs = 7, rowprefix = "", colprefix = "Species")

    @testset "correspondence_analysis" begin
        ca_results = VegSci.correspondence_analysis(x)
        @test typeof(ca_results) <: NamedTuple
        @test keys(ca_results) == (:sv, :rownames, :rowmass, :rowcoord, :colnames, :colmass, :colcoord, :N)
        @test typeof(ca_results.rowcoord) <: NamedMatrix
        @test typeof(ca_results.colcoord) <: NamedMatrix
        @test size(ca_results.rowcoord)[1] == size(x)[1]
        @test size(ca_results.colcoord)[1] == size(x)[2]
        @test size(ca_results.rowcoord)[2] == (size(x)[2] - 1)
        @test size(ca_results.colcoord)[2] == (size(x)[2] - 1)
        @test names(ca_results.rowcoord)[1] == names(x)[1]
        @test names(ca_results.colcoord)[1] == names(x)[2]
    end
end