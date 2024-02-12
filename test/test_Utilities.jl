using EcoVeg
using Test
using NamedArrays

@testset "Utilities.jl" begin
    x = generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "SiteA-", colprefix = "Species")
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