using VegSci
using Documenter
using Test

# Doctests
doctest(VegSci)

# Collate all tests
@testset "VegSci.jl" begin
    include("test_Utilities.jl")
    include("test_SyntopicTables.jl")
    include("test_SimilarityIndices.jl")
    # include("test_CorrespondenceAnalysis.jl")
    # include("test_Fidelity.jl")
    # include("test_Aqua.jl")
end