using VegSci
using Test
using OrderedCollections
using NamedArrays
using SparseArrays
using DataFrames

@testset "Utilities.jl" begin
    test_rowprefix = "Species"
    test_colprefix = "Releve"
    test_rowdim = "Species"
    test_coldim = "Releve"
    @testset "generate_test_array, dense" begin
        test_rown = 10
        test_coln = 10
        test_meancoloccs = 5
        test_val_type = "integer"
        test_scale = "cols"
        test_sparse_array = false
        gtad_results = VegSci.generate_test_array(rown = test_rown, coln = test_coln, 
                                                  meancoloccs = test_meancoloccs, val_type = test_val_type, scale = test_scale,
                                                  rowprefix = test_rowprefix, colprefix = test_colprefix, 
                                                  rowdim = test_rowdim, coldim = test_coldim,
                                                  sparse_array = test_sparse_array)
        @test typeof(gtad_results) <: NamedMatrix{Float64, Matrix{Float64}, Tuple{OrderedDict{String, Int64}, OrderedDict{String, Int64}}}
        @test size(gtad_results) == (test_rown, test_coln)
        @test names(gtad_results)[1] == vec([string("Releve")].*string.([1:1:test_rown;]))
        @test names(gtad_results)[2] == vec([string("Species")].*string.([1:1:test_coln;])) 
        @test all(x-> x >= 0.0, gtad_results)
        @test all(x-> x <= 1.0, gtad_results)
        @test all(sum(gtad_results, dims = 2) .≈ 1.0)
    end
    @testset "generate_test_array, sparse" begin
        test_rown = 1000
        test_coln = 10000
        test_meancoloccs = 15
        test_val_type = "integer"
        test_scale = "cols"
        test_sparse_array = true
        gtas_results = VegSci.generate_test_array(rown = test_rown, coln = test_coln, 
                                                  meancoloccs = test_meancoloccs, val_type = test_val_type, scale = test_scale,
                                                  rowprefix = test_rowprefix, colprefix = test_colprefix, 
                                                  rowdim = test_rowdim, coldim = test_coldim,
                                                  sparse_array = test_sparse_array)
        @test typeof(gtas_results) <: NamedMatrix{Float64, SparseMatrixCSC{Float64, Int64}, Tuple{OrderedDict{String, Int64}, OrderedDict{String, Int64}}}
        @test size(gtas_results) == (test_rown, test_coln)
        @test names(gtas_results)[1] == vec([string("Releve")].*string.([1:1:test_rown;]))
        @test names(gtas_results)[2] == vec([string("Species")].*string.([1:1:test_coln;])) 
        @test all(x-> x >= 0, gtas_results)
        @test all(x-> x <= 100, gtas_results)
        @test all(sum(gtas_results, dims = 2) .>= 90)
        @test all(sum(gtas_results, dims = 2) .<= 110)
    end
    # Generate some arrays to use in the tests below
    x = VegSci.generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "ReleveA-", colprefix = "Species")
    y = VegSci.generate_test_array(rown = 5, coln = 10, meancoloccs = 5, rowprefix = "ReleveB-", colprefix = "Species")
    xs = VegSci.generate_test_array(rown = 100, coln = 100, meancoloccs = 5, rowprefix = "ReleveA-", colprefix = "Species", sparse_array = true)
    ys = VegSci.generate_test_array(rown = 50, coln = 100, meancoloccs = 5, rowprefix = "ReleveB-", colprefix = "Species", sparse_array = true)
    @testset "nm_to_df, dense" begin
        nmtdf_results = VegSci.nm_to_df(x)
        @test typeof(nmtdf_results) <: DataFrame
        @test nmtdf_results[!, 1] == names(x, 1)
        @test names(nmtdf_results) == [dimnames(x, 1); names(x, 2)]
    end
    @testset "nm_to_df, sparse" begin
        nmtdf_results = VegSci.nm_to_df(xs)
        @test typeof(nmtdf_results) <: DataFrame
        @test nmtdf_results[!, 1] == names(xs, 1)
        @test names(nmtdf_results) == [dimnames(xs, 1); names(xs, 2)]
    end
    @testset "nzfunc, dense" begin
        nzfunc_results = VegSci.nzfunc(minimum, x, dims = 1) 
        @test typeof(nzfunc_results) <: NamedMatrix
        @test size(nzfunc_results) == (1, 10)
        @test all(nzfunc_results != 0)
        @test names(nzfunc_results)[2] == names(x)[2]
    end
    @testset "merge_namedarrays, dense" begin
        merna_results = VegSci.merge_namedarrays([x[:,Not(["Species3", "Species10"])], y[:,Not(["Species4"])]])
        @test typeof(merna_results) <: NamedMatrix
        @test size(merna_results) == (15, 10)
        @test names(merna_results)[1] == cat(names(x)[1], names(y)[1], dims = 1) 
        @test Set(names(merna_results)[2]) == Set(unique(cat(names(x)[2], names(y)[2], dims = 1)))
    end
    # @testset "merge_namedarrays, sparse" begin
    #     merna_results = VegSci.merge_namedarrays([xs[:,Not(["Species3", "Species10"])], ys[:,Not(["Species4"])]])
    #     @test typeof(merna_results) <: NamedMatrix
    #     @test size(merna_results) == (15, 10)
    #     @test names(merna_results)[1] == cat(names(x)[1], names(y)[1], dims = 1) 
    #     @test Set(names(merna_results)[2]) == Set(unique(cat(names(x)[2], names(y)[2], dims = 1)))
    # end
    @testset "align_array_columns, dense" begin
        aac_results = VegSci.align_array_columns(x[:,Not(["Species3", "Species10"])], y[:,Not(["Species4"])])
        @test names(aac_results.x)[2] == names(aac_results.y)[2]
    end
    @testset "align_array_columns, sparse" begin
        aac_results = VegSci.align_array_columns(xs[:,Not(["Species3", "Species10"])], ys[:,Not(["Species4"])])
        @test names(aac_results.x)[2] == names(aac_results.y)[2]
    end
end