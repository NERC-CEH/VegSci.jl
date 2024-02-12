using EcoVeg
using Test
using Suppressor
using NamedArrays
using DataFrames

@testset "SyntopicTables.jl" begin
    x = generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "SiteA-", colprefix = "Species")
    csto_results = EcoVeg.compose_syntopic_table_object("Test", x)
    @testset "compose_syntopic_table_object" begin 
        @test typeof(csto_results) <: SyntopicTable
        @test fieldnames(typeof(csto_results)) == (:name, :table, :releve_ids, :releve_n, :species_n, :abundance_units)
        @test typeof(csto_results.name) <: String
        @test typeof(csto_results.table) <: DataFrame
        @test typeof(csto_results.releve_ids) <: Vector{String}
        @test typeof(csto_results.releve_n) <: Int64
        @test typeof(csto_results.species_n) <: Int64
        @test typeof(csto_results.abundance_units) <: String
    end
    @testset "print_summary_syntopic_table" begin
        # Need to think about how I test this given the complexity of the output.
        psst_result = @capture_out EcoVeg.print_summary_syntopic_table(csto_results)
        @test typeof(psst_result) <: String
    end
end