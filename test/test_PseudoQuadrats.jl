using VegSci
using Test
using NamedArrays
using DataFrames

@testset "PseudoQuadrats.jl" begin
    ta = VegSci.generate_test_array(rown = 15, coln = 30, meancoloccs = 10, rowprefix = "SiteA-", colprefix = "Species")
    sto = VegSci.compose_syntopic_table_object("Test", "T", ta)
    n = 15
    @testset "generate_psquads_sto" begin
        gpsqsto_results  = VegSci.generate_psquads_sto(sto, n)
        @test typeof(gpsqsto_results) <: NamedMatrix
        @test size(gpsqsto_results) == (15, 30)
        @test all(x -> x in [0, 1], gpsqsto_results)
        @test names(gpsqsto_results)[1] == string.(1:n)
        @test names(gpsqsto_results)[2] == sto.species_names
    end
    @testset "generate_psquads" begin
        species = sto.species_names
        probs = sto.relative_frequency
        gpsq_results  = VegSci.generate_psquads(species, probs, n)
        @test typeof(gpsq_results) <: NamedMatrix
        @test size(gpsq_results) == (15, 30)
        @test all(x -> x in [0, 1], gpsq_results)
        @test names(gpsq_results)[1] == string.(1:n)
        @test names(gpsq_results)[2] == species
    end
end