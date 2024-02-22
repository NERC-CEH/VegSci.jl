using VegSci
using Test
using NamedArrays

@testset "Utilities.jl" begin
    x = generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "SiteA-", colprefix = "Species")
    y = generate_test_array(rown = 5, coln = 10, meancoloccs = 5, rowprefix = "SiteB-", colprefix = "Species")
    @testset "generate_test_array" begin
        gta_results = VegSci.generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "Releve", colprefix = "Species")
        @test typeof(gta_results) <: NamedMatrix
        @test size(gta_results) == (10, 10)
        @test names(gta_results)[1] == vec([string("Releve")].*string.([1:1:10;]))
        @test names(gta_results)[2] == vec([string("Species")].*string.([1:1:10;])) 
        @test all(x->x>=0.0, gta_results)
        @test all(x->x<=1.0, gta_results)
        @test all(sum(x, dims = 2) .≈ 1.0)
    end
    @testset "nzfunc" begin
        nzfunc_results = VegSci.nzfunc(minimum, x, dims = 1) 
        @test typeof(nzfunc_results) <: NamedMatrix
        @test size(nzfunc_results) == (1, 10)
        @test all(nzfunc_results != 0)
        @test names(nzfunc_results)[2] == names(x)[2]
    end
    @testset "merge_namedarrays" begin
        merna_results = VegSci.merge_namedarrays([x[:,Not(["Species3", "Species10"])], y[:,Not(["Species4"])]])
        @test typeof(merna_results) <: NamedMatrix
        @test size(merna_results) == (15, 10)
        @test names(merna_results)[1] == cat(names(x)[1], names(y)[1], dims = 1) 
        @test Set(names(merna_results)[2]) == Set(unique(cat(names(x)[2], names(y)[2], dims = 1)))
    end
end