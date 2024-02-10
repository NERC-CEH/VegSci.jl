using EcoVeg
using NamedArrays
using Suppressor
using DataFrames
using Test

x = generate_test_array(rown = 10, coln = 10, min = 0.0, max = 1.0, increment = 0.1, rowprefix = "Releve", colprefix = "Species")
y = generate_test_array(rown = 5, coln = 10, min = 0.0, max = 1.0, increment = 0.1, rowprefix = "Releve", colprefix = "Species")


@testset "EcoVeg.jl" begin

    @testset "Utilities.jl" begin
        @testset "generate_test_array" begin
            gta_results  = EcoVeg.generate_test_array(rown = 10, coln = 10, min = 0.0, max = 1.0, increment = 0.1, rowprefix = "Releve", colprefix = "Species")
            @test typeof(gta_results) <: NamedMatrix
            @test size(gta_results) == (10, 10)
            @test names(gta_results)[1] == vec([string("Releve")].*string.([1:1:10;]))
            @test names(gta_results)[2] == vec([string("Species")].*string.([1:1:10;])) 
            @test all(x->x>=0.0, gta_results)
            @test all(x->x<=1.0, gta_results)
        end
        @testset "align_array_columns" begin
            aac_results = EcoVeg.align_array_columns(x[:,Not(["Species3", "Species10"])], y[:,Not(["Species4"])])
            @test keys(aac_results) == (:x, :y)
            @test names(aac_results.x)[2] == names(aac_results.y)[2]
        end
    end
    @testset "SyntopicTables.jl" begin
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
            psst_result = @capture_out print_summary_syntopic_table(csto_results)
        end
    end
    @testset "SimilarityIndices.jl" begin
        @testset begin
            binsim_results = EcoVeg.binary_similarity(x, "(a ./ (a .+ b .+ c)) + I")
        end
        @testset begin
            czeksim_results = EcoVeg.czekanowski_index(x, y)
            @test typeof(czeksim_results) <: NamedMatrix
            @test size(czeksim_results) == (10, 5)
            @test all(x->x>=0.0, czeksim_results)
            @test all(x->x<=1.0, czeksim_results)
            @test names(czeksim_results)[1] == names(x)[1]
            @test names(czeksim_results)[2] == names(y)[1]
        end
    end
end
