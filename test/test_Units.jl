using VegSci
using Test
using NamedArrays

@testset "Units.jl" begin

    x = VegSci.generate_test_array(rown = 30, coln = 10, meancoloccs = 7, rowprefix = "", colprefix = "Species")

    @testset "prop_to_domin" begin
        proptd_results = VegSci.prop_to_domin(x)
        @test typeof(proptd_results) <: NamedMatrix
        @test issubset(Set(unique(proptd_results)), Set([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]))
    end
    @testset "perc_to_domin" begin
        perctd_results = VegSci.perc_to_domin(x)
        @test typeof(perctd_results) <: NamedMatrix
        @test issubset(Set(unique(perctd_results)), Set([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]))
    end
    @testset "relfreq_to_constancy" begin
        rftc_results = VegSci.relfreq_to_constancy(x)
        @test typeof(rftc_results) <: NamedMatrix
        @test issubset(Set(unique(rftc_results)), Set([1.0, 2.0, 3.0, 4.0, 5.0]))
    end
end