using VegSci
using Test
using DataFrames

@testset "Validation.jl" begin

    clusters_t = Dict("1" => ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
                      "2" => ["11", "12", "13", "14", "15", "16", "17", "18", "19", "20"],
                      "3" => ["21", "22", "23", "24", "25", "26", "27", "28", "29", "30"])
    clusters_f = Dict("1" => ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
                      "2" => ["1", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"],
                      "3" => ["17", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"])

    @testset "check_cluster_releveID_duplicates" begin
        ccridd_results_t = VegSci.check_cluster_releveID_duplicates(clusters_t)
        @test typeof(ccridd_results_t) <: DataFrame
        @test nrow(ccridd_results_t) == 0
        ccridd_results_f = VegSci.check_cluster_releveID_duplicates(clusters_f)
        @test typeof(ccridd_results_f) <: DataFrame
        @test nrow(ccridd_results_f) == 2
        @test ccridd_results_f.duplicates == vec(["1", "17"])
    end
end