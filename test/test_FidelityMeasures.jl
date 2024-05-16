using VegSci
using Test
using NamedArrays
using DataFrames

@testset "FidelityMeasures.jl" begin

    x = VegSci.generate_test_array(rown = 30, coln = 10, meancoloccs = 7, rowprefix = "", colprefix = "Species")
    # Set clusters manually for now as I don't yet want to commit to making clustering a dependency.
    clusters = Dict("1" => ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
                    "2" => ["11", "12", "13", "14", "15", "16", "17", "18", "19", "20"],
                    "3" => ["21", "22", "23", "24", "25", "26", "27", "28", "29", "30"])
    csto_results = VegSci.compose_syntopic_table_object("Test", "T", x)

    @testset "indval_fidelity" begin
        indval_results = VegSci.indval_fidelity(x, clusters)
        @test typeof(indval_results) <: NamedMatrix
        @test size(indval_results) == (3, 10)
        @test all(x -> x .>= 0.0, indval_results)
        @test all(x -> x .<= 1.0, indval_results)
        @test names(indval_results)[1] == names(clusters)
        @test names(indval_results)[2] == names(x)[2]
    end
    @testset "u_fidelity" begin
        @testset "u_hyp" begin
            u_hyp_results = VegSci.u_fidelity(x, clusters, "u_hyp")
            @test typeof(u_hyp_results) <: NamedMatrix
            @test size(u_hyp_results) == (3, 10)
            @test names(u_hyp_results)[1] == names(clusters)
            @test names(u_hyp_results)[2] == names(x)[2]
        end
        @testset "u_binA" begin
            u_binA_results = VegSci.u_fidelity(x, clusters, "u_binA")
            @test typeof(u_binA_results) <: NamedMatrix
            @test size(u_binA_results) == (3, 10)
            @test names(u_binA_results)[1] == names(clusters)
            @test names(u_binA_results)[2] == names(x)[2]
        end
        @testset "u_binB" begin
            u_binB_results = VegSci.u_fidelity(x, clusters, "u_binB")
            @test typeof(u_binB_results) <: NamedMatrix
            @test size(u_binB_results) == (3, 10)
            @test names(u_binB_results)[1] == names(clusters)
            @test names(u_binB_results)[2] == names(x)[2]
        end
    end
    @testset "phi_fidelity" begin
        phi_results = VegSci.phi_fidelity(x, clusters)
        @test typeof(phi_results) <: NamedMatrix
        @test size(phi_results) == (3, 10)
        @test all(x -> x .>= -1.0, phi_results)
        @test all(x -> x .<= 1.0, phi_results)
        @test names(phi_results)[1] == names(clusters)
        @test names(phi_results)[2] == names(x)[2]
    end
    @testset "chisq_fidelity" begin
        @testset "chisq" begin
           chisq_results = VegSci.chisq_fidelity(x, clusters, "chisq")
           @test typeof(chisq_results) <: NamedMatrix
           @test size(chisq_results) == (3, 10)
           @test names(chisq_results)[1] == names(clusters)
           @test names(chisq_results)[2] == names(x)[2]
        end
        @testset "chisq_adj" begin
           chisq_adj_results = VegSci.chisq_fidelity(x, clusters, "chisq_adj")
           @test typeof(chisq_adj_results) <: NamedMatrix
           @test size(chisq_adj_results) == (3, 10)
           @test names(chisq_adj_results)[1] == names(clusters)
           @test names(chisq_adj_results)[2] == names(x)[2]
        end
    end
    @testset "G_fidelity" begin
        @testset "G" begin
           G_results = VegSci.G_fidelity(x, clusters, "G")
           @test typeof(G_results) <: NamedMatrix
           @test size(G_results) == (3, 10)
           @test names(G_results)[1] == names(clusters)
           @test names(G_results)[2] == names(x)[2]
        end
        @testset "G_adj" begin
           G_adj_results = VegSci.G_fidelity(x, clusters, "G_adj")
           @test typeof(G_adj_results) <: NamedMatrix
           @test size(G_adj_results) == (3, 10)
           @test names(G_adj_results)[1] == names(clusters)
           @test names(G_adj_results)[2] == names(x)[2]
        end
    end
    @testset "fidelity_utlities" begin
        phi_results = VegSci.phi_fidelity(x, clusters)
        ei_results = VegSci.extract_indicators(phi_results; p_cut = 0.3, n_cut = -0.1, tp = true)
        @testset "extract_indicators" begin
            @test typeof(ei_results) <: NamedMatrix
        end
        @testset "ind_mat_to_df" begin
            imtdf_results = VegSci.ind_mat_to_df(ei_results)
            @test typeof(imtdf_results) <: DataFrame
        end
        @testset "assign_fidelity_synobj" begin
            afs_results = VegSci.assign_fidelity_synobj(csto_results, phi_results, 0.3, -0.1)
            @test typeof(afs_results) <: SyntopicTable
            @test typeof(afs_results.minimum_species) <: Int64
            @test typeof(afs_results.mean_species) <: Float64
            @test typeof(afs_results.maximum_species) <: Int64
        end 
    end
end