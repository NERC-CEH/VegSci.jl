module VegSci

export generate_test_array
export nzfunc
export merge_namedarrays
export align_array_columns
include("Utilities.jl")

export check_cluster_releveID_duplicates
include("Validation.jl")

export SyntopicTable
export compose_syntopic_table_object
export extract_syntopic_matrix
export extract_syntopic_table
export print_summary_syntopic_table
include("SyntopicTables.jl")

export indval_fidelity
export u_fidelity
export phi_fidelity
export chisq_fidelity
export G_fidelity
export extract_indicators
export ind_mat_to_df
export assign_fidelity_synobj
include("FidelityMeasures.jl")

export jaccard_coefficient
export steinhaus_coefficient
# export similarity_fqi
include("SimilarityIndices.jl")

export correspondence_analysis
include("CorrespondenceAnalysis.jl")

end
