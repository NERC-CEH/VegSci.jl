using VegSci
using Test
using DataFrames

@testset "Validation.jl" begin
    x = VegSci.generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "SiteA-", colprefix = "Species")
    x_f = NamedArray(rand(-2.0:0.1:2.0, 5, 5), names = (["Releve1", "Releve2", 3, "Releve4", "Releve5"], ["Spp1", "Spp2", "Spp3", "Spp4", "Spp5"]))
    x_f["Releve1", "Spp4"] = NaN
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
    @testset "check_relSpp_array_format" begin
        crsaf_results_t = VegSci.check_relSpp_array_format(x)
        @test crsaf_results_t["correct_format"] == true
        crsaf_results_f = VegSci.check_relSpp_array_format(x_f)
        @test crsaf_results_f["correct_format"] == false
    end
end