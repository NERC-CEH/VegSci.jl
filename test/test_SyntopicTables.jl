using VegSci
using Test
using Suppressor
using NamedArrays
using DataFrames

@testset "SyntopicTables.jl" begin
    x = VegSci.generate_test_array(rown = 15, coln = 10, meancoloccs = 7, rowprefix = "SiteA-", colprefix = "Species")
    @testset "compose_syntopic_table_object with no cover_abundance_scale and excl_thresh arguments" begin
        csto_results = VegSci.compose_syntopic_table_object("Test", "T", x)
        @test typeof(csto_results) <: SyntopicTable
        @test fieldnames(typeof(csto_results)) == fieldnames(VegSci.SyntopicTable)
        @test csto_results.additional_species == [""]
    end
    csto_results = VegSci.compose_syntopic_table_object("Test", "T", x, excl_thresh = 0.7)
    @testset "compose_syntopic_table_object" begin 
        @test typeof(csto_results) <: SyntopicTable
        @test fieldnames(typeof(csto_results)) == fieldnames(VegSci.SyntopicTable) 
        @test typeof(csto_results.name) <: String
        @test typeof(csto_results.code) <: String
        @test typeof(csto_results.releve_n) <: Int64
        @test typeof(csto_results.releve_ids) <: Vector{String}
        @test typeof(csto_results.species_n) <: Int64
        @test typeof(csto_results.species_names) <: Vector{String}
        @test typeof(csto_results.minimum_species) <: Int64
        @test typeof(csto_results.mean_species) <: Float64
        @test typeof(csto_results.maximum_species) <: Int64
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
        @test typeof(csto_results.additional_species) <: Vector{String}
    end
    @testset "extract_syntopic_matrix" begin
        esm_results = VegSci.extract_syntopic_matrix(csto_results)
        @test typeof(esm_results) <: NamedMatrix
        @test Set(names(esm_results)[2]) == Set(setdiff(names(x)[2], csto_results.additional_species))
        @test size(esm_results) == (1, csto_results.species_n)
        @test all(x -> x .>= 0.0, esm_results)
        @test all(x -> x .<= 1.0, esm_results)
    end
    @testset "extract_syntopic_table" begin
        est_results = VegSci.extract_syntopic_table(csto_results)
        @test typeof(est_results) <: DataFrame
    end
    @testset "print_summary_syntopic_table" begin
        # Need to think about how I test this given the complexity of the output.
        psst_result = @capture_out VegSci.print_summary_syntopic_table(csto_results, "normal", "proportion")
        @test typeof(psst_result) <: String
    end
end