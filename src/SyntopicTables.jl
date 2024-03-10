using DataFrames
using Statistics
using PrettyTables
using NamedArrays

"""
SyntopicTable

# Fields
- name -- The names of the association, of class String
- releve_n -- The number of releves used to compose the syntopic table object, of class Int64
- releve_ids -- The ids of the releves used in the construction of the syntopic table object, usually the row names of the releve by species matrix, of class Vector{String}
- species_n -- The number of species in the syntopic table, of class Int64
- species_names -- The names of species in the syntopic table, of class Vector{String}
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

A object of class VegSci.SyntopicTable. See ...

### Notes


### Algorithm


### References

"""
function compose_syntopic_table_object(name::String, x::NamedMatrix, cover_abundance_scale::String = "proportion")

    x = x[:, vec(map(col -> any(col .!= 0), eachcol(x)))]
    x = x[vec(map(col -> any(col .!= 0), eachrow(x))), :]

    releve_n = size(x)[1]
    releve_ids = names(x)[1]
    species_n = size(x)[2]
    species_names = names(x)[2]

    abs_frequency = vec(sum(x->x>0, x, dims=1))
    rel_frequency = vec(abs_frequency ./ releve_n)
    min_abundance = vec(VegSci.nzfunc(minimum, x, dims = 1))
    max_abundance = vec(maximum(x, dims = 1))
    mean_abundance = vec(mean(x, dims = 1))
    median_abundance = vec(VegSci.nzfunc(median, x, dims = 1))
    fidelity = [0.0]
    fidelity_p = [""]
    fidelity_n = [""]
    
    syntopic_table_object = VegSci.SyntopicTable(name,
                                                 releve_n,
                                                 releve_ids, 
                                                 species_n, 
                                                 species_names,
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

function extract_syntopic_table(syntopic_table_object::SyntopicTable, frequency_scale::String, cover_abundance_scale::String)

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
    sort!(table, [:RelativeFrequency,:MedianAbundance], rev = true)

    # Convert frequency values to the desired format
    if frequency_scale == "normal"
    end

    # Convert cover-abundance values to the desired format
    if cover_abundance_scale == "proportion"
    end
    
    return table

end

function print_summary_syntopic_table(syntopic_table_object::SyntopicTable, frequency_scale::String, cover_abundance_scale::String)

    table = VegSci.extract_syntopic_table(syntopic_table_object, frequency_scale, cover_abundance_scale)

    summary_syntopic_table = table[:, [:Species, :RelativeFrequency, :Abundance]]

    name = syntopic_table_object.name
    releve_n = syntopic_table_object.releve_n
    species_n = syntopic_table_object.species_n
    fidelity_p = syntopic_table_object.fidelity_p
    fidelity_n = syntopic_table_object.fidelity_n

    println("\n")
    println("Community Name: $name")
    println("Releves: n = $releve_n")
    println("Species: n = $species_n")
    printstyled(string("Postive Indicators: ", join(fidelity_p, " "), "\n"); color = :green)
    printstyled(string("Negative Indicators: ", join(fidelity_n, " "), "\n"); color = :red)
    pretty_table(summary_syntopic_table, show_subheader = false)

end