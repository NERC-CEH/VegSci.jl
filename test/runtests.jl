using EcoVeg
using NamedArrays
using Suppressor
using DataFrames
using Test
using Aqua

@testset "EcoVeg.jl" begin

    x = generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "SiteA-", colprefix = "Species")
    y = generate_test_array(rown = 5, coln = 10, meancoloccs = 5, rowprefix = "SiteB-", colprefix = "Species")

    @testset "Utilities.jl" begin
        @testset "generate_test_array" begin
            gta_results = EcoVeg.generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "Releve", colprefix = "Species")
            @test typeof(gta_results) <: NamedMatrix
            @test size(gta_results) == (10, 10)
            @test names(gta_results)[1] == vec([string("Releve")].*string.([1:1:10;]))
            @test names(gta_results)[2] == vec([string("Species")].*string.([1:1:10;])) 
            @test all(x->x>=0.0, gta_results)
            @test all(x->x<=1.0, gta_results)
            @test all(sum(x, dims = 2) .≈ 1.0)
        end
        @testset "nzfunc" begin
            nzfunc_results = EcoVeg.nzfunc(minimum, x, dims = 1) 
            @test typeof(nzfunc_results) <: NamedMatrix
            @test size(nzfunc_results) == (1, 10)
            @test all(nzfunc_results != 0)
            @test names(nzfunc_results)[2] == names(x)[2]
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
            psst_result = @capture_out EcoVeg.print_summary_syntopic_table(csto_results)
            @test typeof(psst_result) <: String
        end
    end
    @testset "SimilarityIndices.jl" begin
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

# @testset "Aqua.jl" begin
#     Aqua.test_all(
#       EcoVeg;
#       piracies=false,
#     )
#   end