using VegSci
using Test
using Suppressor
using NamedArrays
using DataFrames

@testset "SyntopicTables.jl" begin
    x = VegSci.generate_test_array(rown = 10, coln = 10, meancoloccs = 7, rowprefix = "SiteA-", colprefix = "Species")
    csto_results = VegSci.compose_syntopic_table_object("Test", x)
    @testset "compose_syntopic_table_object" begin 
        @test typeof(csto_results) <: SyntopicTable
        @test fieldnames(typeof(csto_results)) == fieldnames(VegSci.SyntopicTable) 
        @test typeof(csto_results.name) <: String
        @test typeof(csto_results.releve_n) <: Int64
        @test typeof(csto_results.releve_ids) <: Vector{String}
        @test typeof(csto_results.species_n) <: Int64
        @test typeof(csto_results.species_names) <: Vector{String}
        @test typeof(csto_results.abundance_units) <: String
        @test typeof(csto_results.relative_frequency) <: Vector{Float64}
        @test typeof(csto_results.absolute_frequency) <: Vector{Float64}
        @test typeof(csto_results.minimum_abundance) <: Vector{Float64}
        @test typeof(csto_results.mean_abundance) <: Vector{Float64}
        @test typeof(csto_results.median_abundance) <: Vector{Float64}
        @test typeof(csto_results.maximum_abundance) <: Vector{Float64}
        @test typeof(csto_results.fidelity) <: Vector{Float64}
        @test typeof(csto_results.fidelity_p) <: Vector{String}
        @test typeof(csto_results.fidelity_n) <: Vector{String}
    end
    @testset "extract_syntopic_matrix" begin
        esm_results = VegSci.extract_syntopic_matrix(csto_results)
        @test typeof(esm_results) <: NamedMatrix
        @test names(esm_results)[2] == names(x)[2]
        @test size(esm_results) == (1, 10)
        @test all(x -> x .>= 0.0, esm_results)
        @test all(x -> x .<= 1.0, esm_results)
    end
    @testset "extract_syntopic_table" begin
        est_results = VegSci.extract_syntopic_table(csto_results)
        @test typeof(est_results) <: DataFrame
    end
    @testset "print_summary_syntopic_table" begin
        # Need to think about how I test this given the complexity of the output.
        psst_result = @capture_out VegSci.print_summary_syntopic_table(csto_results, "normal", "cover_proportion")
        @test typeof(psst_result) <: String
    end
end