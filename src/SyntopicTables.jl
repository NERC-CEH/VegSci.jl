using DataFrames
using Statistics
using PrettyTables

struct SyntopicTable
    name::String
    table::DataFrame
    releve_ids::Vector{String}
    releve_n::Int64
    species_n::Int64
    abundance_units::String
end

function compose_syntopic_table_object(name::String, mat::NamedArray)

    releve_n = size(mat)[1]
    releve_ids = names(mat)[1]
    species_names = names(mat)[2]
    species_n = size(mat)[2]
    abs_frequency = sum(x->x>0, mat, dims=1)
    rel_frequency = abs_frequency ./ releve_n
    mat_nozeros = copy(mat)
    mat_nozeros[iszero.(mat_nozeros)] .= 1
    min_abundance = minimum(mat_nozeros, dims = 1)
    max_abundance = maximum(mat, dims = 1)
    mean_abundance = mean(mat, dims = 1)
    median_abundance = median(mat, dims = 1)

    # Create table
    table = DataFrame(Species = vec(species_names), RelativeFrequency = vec(rel_frequency), AbsoluteFrequency = vec(abs_frequency),
                      MinimumAbundance = vec(min_abundance), MaximumAbundance = vec(max_abundance), MeanAbundance = vec(mean_abundance),
                      MedianAbundance = vec(median_abundance))
    table[!, :Abundance] = string.(table.MedianAbundance, " (", table.MinimumAbundance, " - ", table.MaximumAbundance, ")")

    # Order table by relative frequency
    sort!(table, [:RelativeFrequency, :MedianAbundance], rev = [true, true])
    
    # Form object
    syntopic_table_object = SyntopicTable(name, table, releve_ids, releve_n, species_n, "cover_proportion")

    return syntopic_table_object

end

function print_summary_syntopic_table(syntopic_table_object::SyntopicTable, frequency_scale::String = "normal", cover_abundance_scale::String = "cover_proportion")

    # Select relevant columns
    summary_syntopic_table = syntopic_table_object.table[:, [:Species, :RelativeFrequency, :Abundance]]

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