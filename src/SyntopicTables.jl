using DataFrames
using Statistics
using PrettyTables
using NamedArrays

"""
SyntopicTable

# Fields
- name -- The name of the association, of class String
- code -- The code for the association, of class String
- releve_n -- The number of releves used to compose the syntopic table object, of class Int64
- releve_ids -- The ids of the releves used in the construction of the syntopic table object, usually the row names of the releve by species matrix, of class Vector{String}
- species_n -- The number of species in the syntopic table, of class Int64
- species_names -- The names of species in the syntopic table, of class Vector{String}
- minimum_species -- The minimum species richness within the releves, of class Int64
- mean_species -- The mean species richness across the releves, of class Float64
- maximum_species -- The maximum species richness within the releves, of class Int64
- abundance_units -- The abundance units, of class String
- relative_frequency -- The relative frequency of occurrence of each species across all releves, of class Vector{Float64}
- absolute_frequency -- Vector{Float64}
- minimum_abundance -- Vector{Float64}
- mean_abundance -- Vector{Float64}
- median_abundance -- Vector{Float64}
- maximum_abundance -- Vector{Float64}
- fidelity -- Vector{Float64}
- fidelity_p -- Vector{String}
- fidelity_n -- Vector{String}
"""
mutable struct SyntopicTable
    name::String
    code::String
    releve_n::Int64
    releve_ids::Vector{String}
    species_n::Int64
    species_names::Vector{String}
    minimum_species::Int64
    mean_species::Float64
    maximum_species::Int64
    abundance_units::String
    relative_frequency::Vector{Float64}
    absolute_frequency::Vector{Float64}
    minimum_abundance::Vector{Float64}
    mean_abundance::Vector{Float64}
    median_abundance::Vector{Float64}
    maximum_abundance::Vector{Float64}
    fidelity::Vector{Float64}
    fidelity_p::Vector{String}
    fidelity_n::Vector{String}
end

"""
compose_syntopic_table_object(name::String, x::NamedArray)::SyntopicTable

Compose a single syntopic table object from a releve by species matrix 

### Input

- `name` -- The name given to the association represented in the sytopic table, of class String
- `x` -- A releve by species matrix, of class NamedArrays::NamedMatrix
- `cover_abundance_scale` -- The units of the matrix values, of class String

### Output

A object of class VegSci.SyntopicTable.

### Notes


### Algorithm


### References

"""
function compose_syntopic_table_object(name::String, code::String, x::NamedMatrix, cover_abundance_scale::String = "proportion")

    # Check the dimension names of the input matrix are vectors of strings
    @assert typeof(names(x)) <: Vector{Vector{String}} "Both the row names and column names must only contain strings."

    # Remove columns containing all zeros (species which do not occur in any releves)
    x = x[:, vec(map(col -> any(col .!= 0), eachcol(x)))]
    x = x[vec(map(col -> any(col .!= 0), eachrow(x))), :]

    species_count = vec(sum(x->x>0, x, dims=2))

    releve_n = size(x)[1]
    releve_ids = names(x)[1]
    species_n = size(x)[2]
    species_names = names(x)[2]
    minimum_species = minimum(species_count)
    mean_species = mean(species_count)
    maximum_species = maximum(species_count)
    abs_frequency = vec(sum(x->x>0, x, dims=1))
    rel_frequency = vec(abs_frequency ./ releve_n)
    min_abundance = vec(VegSci.nzfunc(minimum, x, dims = 1))
    max_abundance = vec(maximum(x, dims = 1))
    mean_abundance = vec(mean(x, dims = 1))
    median_abundance = vec(VegSci.nzfunc(median, x, dims = 1))
    fidelity = [0.0]
    fidelity_p = [""]
    fidelity_n = [""]
    
    syntopic_table_object = VegSci.SyntopicTable(code,
                                                 name,
                                                 releve_n,
                                                 releve_ids, 
                                                 species_n, 
                                                 species_names,
                                                 minimum_species,
                                                 mean_species,
                                                 maximum_species,
                                                 cover_abundance_scale,
                                                 rel_frequency,
                                                 abs_frequency,
                                                 min_abundance,
                                                 mean_abundance,
                                                 median_abundance,
                                                 max_abundance,
                                                 fidelity,
                                                 fidelity_p,
                                                 fidelity_n
                                                 )

    return syntopic_table_object

end

function extract_syntopic_matrix(syntopic_table_object::SyntopicTable)

    rowname = [syntopic_table_object.name]
    species_names = syntopic_table_object.species_names
    median_abundance = syntopic_table_object.median_abundance
    syntopic_matrix = NamedArrays.NamedArray(reshape(median_abundance, 1, :), names = (rowname, species_names))

    return syntopic_matrix

end

function extract_syntopic_table(syntopic_table_object::SyntopicTable; include_code = false)

    # Create table
    table = DataFrame(Species = syntopic_table_object.species_names, 
                      RelativeFrequency = syntopic_table_object.relative_frequency, 
                      AbsoluteFrequency = syntopic_table_object.absolute_frequency,
                      MinimumAbundance = syntopic_table_object.minimum_abundance, 
                      MaximumAbundance = syntopic_table_object.maximum_abundance, 
                      MeanAbundance = syntopic_table_object.mean_abundance,
                      MedianAbundance = syntopic_table_object.median_abundance
                      )

    if include_code == true
        insertcols!(table, 1, :CommunityCode => syntopic_table_object.code)
    end

    # Order table by relative frequency
    sort!(table, [:RelativeFrequency,:MedianAbundance], rev = true)
    
    return table

end

function print_summary_syntopic_table(syntopic_table_object::SyntopicTable, frequency_scale::String, cover_abundance_scale::String)

    # Convert frequency values to the desired format
    if frequency_scale == "normal"
    end

    # Convert cover-abundance values to the desired format
    if cover_abundance_scale == "proportion"
    end

    table = VegSci.extract_syntopic_table(syntopic_table_object)
    table[!, [:RelativeFrequency, :MinimumAbundance, :MeanAbundance, :MaximumAbundance]] = round.(table[:, [:RelativeFrequency, :MinimumAbundance, :MeanAbundance, :MaximumAbundance]], digits = 1)
    table[!, :Abundance] = string.(table.MeanAbundance, " (", table.MinimumAbundance, " - ", table.MaximumAbundance, ")")

    summary_syntopic_table = table[:, [:Species, :RelativeFrequency, :Abundance]]

    name = syntopic_table_object.name
    releve_n = syntopic_table_object.releve_n
    species_n = syntopic_table_object.species_n
    species_min = syntopic_table_object.minimum_species
    species_mean = Int.(round(syntopic_table_object.mean_species, digits = 0))
    species_max = syntopic_table_object.maximum_species
    fidelity_p = syntopic_table_object.fidelity_p
    fidelity_n = syntopic_table_object.fidelity_n

    println("\n")
    println("Community Name: $name")
    println("Releves: n = $releve_n")
    println("Species: n = $species_n")
    println("Species Richness: $species_mean ($species_min - $species_max)")
    printstyled(string("Postive Indicators: ", join(fidelity_p, " "), "\n"); color = :green)
    printstyled(string("Negative Indicators: ", join(fidelity_n, " "), "\n"); color = :red)
    pretty_table(summary_syntopic_table, show_subheader = false)

end