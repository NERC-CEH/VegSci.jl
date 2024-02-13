using DataFrames
using Statistics
using PrettyTables
using NamedArrays

struct SyntopicTable
    name::String
    releve_n::Int64
    releve_ids::Vector{String}
    species_n::Int64
    species_names::Vector{String}
    abundance_units::String
    relative_frequency::Vector{Float64}
    absolute_frequency::Vector{Float64}
    minimum_abundance::Vector{Float64}
    mean_abundance::Vector{Float64}
    median_abundance::Vector{Float64}
    maximum_abundance::Vector{Float64}
end

function compose_syntopic_table_object(name::String, mat::NamedArray)

    mat = mat[:, vec(map(col -> any(col .!= 0), eachcol(mat)))]
    mat = mat[vec(map(col -> any(col .!= 0), eachrow(mat))), :]

    releve_n = size(mat)[1]
    releve_ids = names(mat)[1]
    species_n = size(mat)[2]
    species_names = names(mat)[2]

    abs_frequency = vec(sum(x->x>0, mat, dims=1))
    rel_frequency = vec(abs_frequency ./ releve_n)
    min_abundance = vec(EcoVeg.nzfunc(minimum, mat, dims = 1))
    max_abundance = vec(maximum(mat, dims = 1))
    mean_abundance = vec(mean(mat, dims = 1))
    median_abundance = vec(EcoVeg.nzfunc(median, mat, dims = 1))
    
    syntopic_table_object = SyntopicTable(name,
                                          releve_n,
                                          releve_ids, 
                                          species_n, 
                                          species_names,
                                          "cover_proportion",
                                          rel_frequency,
                                          abs_frequency,
                                          min_abundance,
                                          mean_abundance,
                                          median_abundance,
                                          max_abundance
                                          )

    return syntopic_table_object

end

function print_summary_syntopic_table(syntopic_table_object::SyntopicTable, frequency_scale::String, cover_abundance_scale::String)
        
    # Create table
    table = DataFrame(Species = syntopic_table_object.species_names, 
                      RelativeFrequency = round.(syntopic_table_object.relative_frequency, digits = 1), 
                      AbsoluteFrequency = round.(syntopic_table_object.absolute_frequency, digits = 1),
                      MinimumAbundance = round.(syntopic_table_object.minimum_abundance, digits = 1), 
                      MaximumAbundance = round.(syntopic_table_object.maximum_abundance, digits = 1), 
                      MeanAbundance = round.(syntopic_table_object.mean_abundance, digits = 1),
                      MedianAbundance = round.(syntopic_table_object.median_abundance, digits = 1)
                      )

    table[!, :Abundance] = string.(table.MedianAbundance, " (", table.MinimumAbundance, " - ", table.MaximumAbundance, ")")

    # Order table by relative frequency
    sort!(table,[:RelativeFrequency,:MedianAbundance], rev = true)

    # Select relevant columns
    summary_syntopic_table = table[:, [:Species, :RelativeFrequency, :Abundance]]

    # Convert frequency values to the desired format
    if frequency_scale == "normal"
    end

    # Convert cover-abundance values to the desired format
    if cover_abundance_scale == "cover_proportion"
    end


    name = syntopic_table_object.name
    releve_n = syntopic_table_object.releve_n
    species_n = syntopic_table_object.species_n

    println("\n")
    println("Community Name: $name")
    println("Releves: n = $releve_n")
    println("Species: n = $species_n")
    pretty_table(summary_syntopic_table, show_subheader = false)

end

function extract_syntopic_matrix(syntopic_table_object::SyntopicTable)

    rowname = [syntopic_table_object.name]
    species_names = syntopic_table_object.species_names
    median_abundance = syntopic_table_object.median_abundance
    syntopic_matrix = NamedArrays.NamedArray(reshape(median_abundance, 1, :), names = (rowname, species_names))

    return syntopic_matrix

end