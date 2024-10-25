using VegSci
using Documenter
using Test

# Doctests
doctest(VegSci)

# Collate all tests
@testset "VegSci.jl" begin
    include("test_Units.jl")
    include("test_Utilities.jl")
    include("test_Validation.jl")
    include("test_SyntopicTables.jl")
    include("test_FidelityMeasures.jl")
    include("test_SimilarityIndices.jl")
    include("test_CorrespondenceAnalysis.jl")
    include("test_PseudoQuadrats.jl")
    # include("test_Aqua.jl")
end