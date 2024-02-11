module EcoVeg

export generate_test_array
export nzfunc
export align_array_columns
include("Utilities.jl")

export binary_similarity
export czekanowski_index
include("SimilarityIndices.jl")

# include(Fidelity.jl)

# export correspondence_analysis
# include("CorrespondenceAnalysis.jl")

export SyntopicTable
export compose_syntopic_table_object
export print_summary_syntopic_table
export extract_syntopic_matrix
include("SyntopicTables.jl")

end
