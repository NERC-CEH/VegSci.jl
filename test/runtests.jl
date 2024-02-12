using EcoVeg
using Documenter
using Test

# Doctests
doctest(EcoVeg)

# Collate all tests
@testset "EcoVeg.jl" begin
    include("test_Utilities.jl")
    include("test_SyntopicTables.jl")
    include("test_SimilarityIndices.jl")
    # include("test_CorrespondenceAnalysis.jl")
    # include("test_Fidelity.jl")
    # include("test_Aqua.jl")
end