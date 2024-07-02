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
- additional_species -- Vector{String}
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
    additional_species::Vector{String}
end

"""
compose_syntopic_table_object(name::String, x::NamedArray)::SyntopicTable

Compose a single syntopic table object from a releve by species matrix 

### Input

- `name` -- The name given to the association represented in the sytopic table, of class String
- `x` -- A releve by species matrix, of class NamedArrays::NamedMatrix
- `cover_abundance_scale` -- The units of the matrix values, of class String
- `excl_thresh` -- An optional value representing the lower relative frequency threshold, greater than 0.0 and less than 1.0.
                   Species occurring at a relative frequency of less than this threshold are excluded from the syntopic table and stored
                   in the additional_species field of the VegSci.SyntopicTable object.

### Output

A object of class VegSci.SyntopicTable.

### Notes


### Algorithm


### References

"""
function compose_syntopic_table_object(name::String, code::String, x::NamedMatrix; cover_abundance_scale::String = "proportion", excl_thresh::Union{Nothing, Float64} = nothing)

    # Check that excl_thresh is nothing, or a value greater than 0.0 and less than 1.0.
    @assert (isnothing(excl_thresh) || (excl_thresh > 0.0 && excl_thresh < 1.0)) "If supplied, excl_thresh must be greater than 0.0 and less than 1.0"

    # Check the dimension names of the input matrix are vectors of strings
    @assert typeof(names(x)) <: Vector{Vector{String}} "Both the row names and column names must only contain strings."

    # Check that all values in x are positive
    @assert all(>=(0), x) "All values in the supplied matrix must be greater than or equal to 0."

    # Remove columns containing all zeros (species which do not occur in any releves)
    x = x[:, vec(map(col -> any(col .!= 0), eachcol(x)))]
    x = x[vec(map(col -> any(col .!= 0), eachrow(x))), :]

    # Calculate the number of releves
    releve_n = size(x)[1]

    # Store the names of the releves
    releve_ids = names(x)[1]

    # Calculate the initial absolute frequency of occurrence of each species across all releves
    init_abs_frequency = vec(sum(x -> x > 0, x, dims = 1))

    # Calculate the initial relative frequency of occurrence of each species across all releves
    init_rel_frequency = vec(init_abs_frequency ./ releve_n)

    # If excl_thresh is not nothing
    if !isnothing(excl_thresh)

        # Filter rel_frequency to retrieve the species that occur at a relative frequency over and under excl_thresh
        spp_over_thresh = names(x, 2)[findall(init_rel_frequency .>= excl_thresh)]
        spp_under_thresh = names(x, 2)[findall(init_rel_frequency .< excl_thresh)]

        # Remove the species occuring at a relative frequency of less than excl_thresh from x
        x = x[:, spp_over_thresh]

        # Store the excluded species names
        additional_species = spp_under_thresh
        
    elseif isnothing(excl_thresh)

        additional_species = [""]

    end

    # Re-calculate the absolute frequency of occurrence of each species across all releves
    abs_frequency = vec(sum(x -> x > 0, x, dims = 1))

    # Re-calculate the relative frequency of occurrence of each species across all releves
    rel_frequency = vec(abs_frequency ./ releve_n)

    # Sum the total number of species in each releve - does this occur before or after removing rare species?
    species_count = vec(sum(x -> x > 0, x, dims = 2))

    # Calculate association metrics
    species_n = size(x)[2]
    species_names = names(x)[2]
    minimum_species = minimum(species_count)
    mean_species = mean(species_count)
    maximum_species = maximum(species_count)
    min_abundance = vec(VegSci.nzfunc(minimum, x, dims = 1))
    max_abundance = vec(maximum(x, dims = 1))
    mean_abundance = vec(mean(x, dims = 1))
    median_abundance = vec(VegSci.nzfunc(median, x, dims = 1))
    fidelity = [0.0]
    fidelity_p = [""]
    fidelity_n = [""]

    # Compose the sytnopic table object    
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
                                                 fidelity_n,
                                                 additional_species
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

    x = VegSci.generate_test_array(rown = 15, coln = 10, meancoloccs = 7, rowprefix = "SiteA-", colprefix = "Species")
    syntopic_table_object = VegSci.compose_syntopic_table_object("Test", "T", x, excl_thresh = 0.7)

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