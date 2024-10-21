using VegSci
using Test
using NamedArrays
using SparseArrays

@testset "SimilarityIndices.jl" begin
    x = VegSci.generate_test_array(rown = 10, coln = 10, meancoloccs = 6, rowprefix = "SiteA-", colprefix = "Species", sparse = false)
    y = VegSci.generate_test_array(rown = 5, coln = 10, meancoloccs = 6, rowprefix = "SiteB-", colprefix = "Species", sparse = false)
    xy = VegSci.merge_namedarrays([x, y])
    x_sparse = VegSci.generate_test_array(rown = 10, coln = 10, meancoloccs = 6, rowprefix = "SiteA-", colprefix = "Species", sparse = true)
    y_sparse = VegSci.generate_test_array(rown = 5, coln = 10, meancoloccs = 6, rowprefix = "SiteB-", colprefix = "Species", sparse = true)
    xy_sparse = VegSci.merge_namedarrays([x, y])
    @testset "steinhaus_coefficient" begin
        steinhaus_results = VegSci.steinhaus_coefficient(x, y)
        @test typeof(steinhaus_results) <: NamedMatrix
        @test size(steinhaus_results) == (10, 5)
        @test all(x -> x .>= 0.0, steinhaus_results)
        @test all(x -> x .<= 1.0, steinhaus_results)
        @test names(steinhaus_results)[1] == names(x)[1]
        @test names(steinhaus_results)[2] == names(y)[1]
    end
    @testset "jaccard_coefficient" begin
        jaccard_results = VegSci.jaccard_coefficient(xy)
        @test typeof(jaccard_results) <: NamedMatrix
        @test size(jaccard_results) == (15, 15)
        @test all(x -> x .>= 0.0, jaccard_results)
        @test all(x -> x .<= 1.0, jaccard_results)
        @test names(jaccard_results)[1] == [names(x)[1]; names(y)[1]]
        @test names(jaccard_results)[2] == [names(x)[1]; names(y)[1]]
    end
    # @testset "similarity_tichy" begin
    #     x = VegSci.generate_test_array(rown = 10, coln = 20, meancoloccs = 15, rowprefix = "SiteX-", colprefix = "Species")[:,Not(["Species5", "Species6", "Species20"])]
    #     y = VegSci.generate_test_array(rown = 5, coln = 20, meancoloccs = 15, rowprefix = "SiteY-", colprefix = "Species")[:,Not(["Species3", "Species11", "Species17"])]
    #     z = VegSci.generate_test_array(rown = 5, coln = 20, meancoloccs = 15, rowprefix = "SiteZ-", colprefix = "Species")[:,Not(["Species3", "Species11", "Species17"])]
    #     yz = vcat(y, z)
    #     setnames!(yz, [names(y)[1]; names(z)[1]], 1)

    #     syn_x = VegSci.compose_syntopic_table_object("x", x)
    #     syn_y = VegSci.compose_syntopic_table_object("y", y)
    #     syn_z = VegSci.compose_syntopic_table_object("z", z)
    #     syn_x_mat = VegSci.extract_syntopic_matrix(syn_x)
    #     syn_y_mat = VegSci.extract_syntopic_matrix(syn_y)
    #     syn_z_mat = VegSci.extract_syntopic_matrix(syn_z)
    #     syn_yz_mat = vcat(syn_y_mat, syn_z_mat)
    #     setnames!(syn_yz_mat, [names(syn_y_mat)[1]; names(syn_z_mat)[1]], 1)

    #     clusters = Dict("Y" => names(y)[1], "Z" => names(z)[1])
    #     syn_yz_fidelity = VegSci.u_fidelity(yz, clusters, "u_hyp")

    #     @testset "fqi" begin
    #         fqi_results = VegSci.similarity_fqi(x, y)
    #         @test typeof(fqi_results) <: NamedMatrix
    #         @test all(x -> x .>= 0.0, fqi_results)
    #         @test all(x -> x .<= 100.0, fqi_results)
    #     end
    #     @testset "pfdi" begin
    #         pfdi_results = VegSci.similarity_pfdi(x, y)
    #         @test typeof(pfdi_results) <: NamedMatrix
    #         @test all(x -> x .>= 0.0, pfdi_results)
    #         @test all(x -> x .<= 100.0, pfdi_results)
    #     end
    #     @testset "nfdi" begin
    #         nfdi_results = VegSci.similarity_pfdi(x, y)
    #         @test typeof(nfdi_results) <: NamedMatrix
    #         @test all(x -> x .>= -100.0, nfdi_results)
    #         @test all(x -> x .<= 0.0, nfdi_results)
    #     end
    #     @testset "fpfi" begin
    #         fpfi_results = VegSci.similarity_pfdi(x, y)
    #         @test typeof(fpfi_results) <: NamedMatrix
    #         @test all(x -> x .>= 0.0, fpfi_results)
    #         @test all(x -> x .<= 100.0, fpfi_results)
    #     end
    #     @testset "fgfi" begin
    #         fgfi_results = VegSci.similarity_pfdi(x, y)
    #         @test typeof(fgfi_results) <: NamedMatrix
    #         @test all(x -> x .>= -50.0, fgfi_results)
    #         @test all(x -> x .<= 100.0, fgfi_results)
    #     end    
    # end
end