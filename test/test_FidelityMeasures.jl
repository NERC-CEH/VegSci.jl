using VegSci
using Test
using NamedArrays

@testset "FidelityMeasures.jl" begin

    x = VegSci.generate_test_array(rown = 30, coln = 10, meancoloccs = 7, rowprefix = "", colprefix = "Species")
    # Set clusters manually for now as I don't yet want to commit to making clustering a dependency.
    clusters = Dict("1" => ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
                    "2" => ["11", "12", "13", "14", "15", "16", "17", "18", "19", "20"],
                    "3" => ["21", "22", "23", "24", "25", "26", "27", "28", "29", "30"])

    @testset "indval_fidelity" begin
        indval_results = VegSci.indval_fidelity(x, clusters)
        @test typeof(indval_results) <: NamedMatrix
        @test size(indval_results) == (3, 10)
        @test all(x -> x .>= 0.0, indval_results)
        @test all(x -> x .<= 1.0, indval_results)
        @test names(indval_results)[1] == names(clusters)
        @test names(indval_results)[2] == names(x)[2]
    end
end