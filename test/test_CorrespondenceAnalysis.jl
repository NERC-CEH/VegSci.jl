using VegSci
using Test
using NamedArrays

@testset "CorrespondenceAnalysis.jl" begin
    test_n_axes = 3
    x = VegSci.generate_test_array(rown = 30, coln = 10, meancoloccs = 7, rowprefix = "", colprefix = "Species")
    xs = VegSci.generate_test_array(rown = 1000, coln = 1000, meancoloccs = 15, rowprefix = "", colprefix = "Species", sparse_array = true)
    @testset "correspondence_analysis, dense, 30x10" begin
        @benchmark ca_results = VegSci.correspondence_analysis(x, axes_n = test_n_axes) samples = 10
        @test typeof(ca_results) <: Dict
        @test names(ca_results) == ["StandardCoordsCols", "PrincipleCoordsCols", "StandardCoordsRows", "PrincipleCoordsRows"]
        @test typeof(ca_results["StandardCoordsCols"]) <: NamedMatrix
        @test typeof(ca_results["PrincipleCoordsCols"]) <: NamedMatrix
        @test typeof(ca_results["StandardCoordsRows"]) <: NamedMatrix
        @test typeof(ca_results["PrincipleCoordsRows"]) <: NamedMatrix
        @test size(ca_results["StandardCoordsCols"])[1] == size(x)[2]
        @test size(ca_results["PrincipleCoordsCols"])[1] == size(x)[2]
        @test size(ca_results["StandardCoordsRows"])[1] == size(x)[1]
        @test size(ca_results["PrincipleCoordsRows"])[1] == size(x)[1]
        @test size(ca_results["StandardCoordsCols"])[2] == test_n_axes
        @test size(ca_results["PrincipleCoordsCols"])[2] == test_n_axes
        @test size(ca_results["StandardCoordsRows"])[2] == test_n_axes
        @test size(ca_results["PrincipleCoordsRows"])[2] == test_n_axes
        @test names(ca_results["StandardCoordsCols"])[1] == names(x)
        @test names(ca_results["PrincipleCoordsCols"])[1] == names(x)
        @test names(ca_results["StandardCoordsRows"])[1] == names(x)
        @test names(ca_results["PrincipleCoordsRows"])[1] == names(x)
    end
    @testset "correspondence_analysis, sparse, 30x10" begin
        ca_results = VegSci.correspondence_analysis(xs)
        @test typeof(ca_results) <: Dict
        @test names(ca_results) == ["StandardCoordsCols", "PrincipleCoordsCols", "StandardCoordsRows", "PrincipleCoordsRows"]
        @test typeof(ca_results["StandardCoordsCols"]) <: NamedMatrix
        @test typeof(ca_results["PrincipleCoordsCols"]) <: NamedMatrix
        @test typeof(ca_results["StandardCoordsRows"]) <: NamedMatrix
        @test typeof(ca_results["PrincipleCoordsRows"]) <: NamedMatrix
        @test size(ca_results["StandardCoordsCols"])[1] == size(x)[2]
        @test size(ca_results["PrincipleCoordsCols"])[1] == size(x)[2]
        @test size(ca_results["StandardCoordsRows"])[1] == size(x)[1]
        @test size(ca_results["PrincipleCoordsRows"])[1] == size(x)[1]
        @test size(ca_results["StandardCoordsCols"])[2] == test_n_axes
        @test size(ca_results["PrincipleCoordsCols"])[2] == test_n_axes
        @test size(ca_results["StandardCoordsRows"])[2] == test_n_axes
        @test size(ca_results["PrincipleCoordsRows"])[2] == test_n_axes
        @test names(ca_results["StandardCoordsCols"])[1] == names(x)
        @test names(ca_results["PrincipleCoordsCols"])[1] == names(x)
        @test names(ca_results["StandardCoordsRows"])[1] == names(x)
        @test names(ca_results["PrincipleCoordsRows"])[1] == names(x)
    end
end