module VegSci

export generate_test_array
export nzfunc
export merge_namedarrays
include("Utilities.jl")

export jaccard_coefficient
export steinhaus_coefficient
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