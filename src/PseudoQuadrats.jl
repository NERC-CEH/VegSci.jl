using Distributions
using NamedArrays

"""
     generate_psquads_sto(sto:SyntopicTable, n::Int64)::NamedMatrix

Generate a set of n pseudo-quadrats from a syntopic table object.

### Input

- `sto` -- A syntopic table object, of class VegSci.SyntopicTable, produced by the `VegSci.compose_syntopic_table_object` function.
- `n` -- The number of pseudo-quadrats to generate, of class Int64.

### Output

A 2 dimensional named array, of class NamedArrays.NamedMatrix, containing pseudo-quadrats for the supplied associations syntopic table, with 1's represent presences, and 0's representing absences.

### Notes

### Algorithm

This function performs a Bernoulli draw using the relative frequencies for each species present in the syntopic table object n times.

### References


"""
function generate_psquads_sto(sto::SyntopicTable, n::Int64)

    # Retrieve the species names and relative frequencies from the syntopic tables object
    species = sto.species_names
    probs = sto.relative_frequency

    # Prepare empty matrix
    psquads = NamedArrays.NamedArray(fill(0, n, length(species)), names = (string.(1:n), species))

    # Create a pseudo-quadrat n times and populate the matrix with the presences
    for i in 1:n
        trial = reduce(vcat, map(x -> rand(Distributions.Bernoulli(x), 1), probs))
        selected_species = species[trial]
        psquads[i, selected_species] .= 1
    end

    return psquads

end

function generate_psquads(species::Vector{String}, probs::Vector{Float64}, n::Int64)

    # Prepare empty matrix
    psquads = NamedArrays.NamedArray(fill(0, n, length(species)), names = (string.(1:n), species))

    # Create a pseudo-quadrat n times and populate the matrix with the presences
    for i in 1:n
        trial = reduce(vcat, map(x -> rand(Distributions.Bernoulli(x), 1), probs))
        selected_species = species[trial]
        psquads[i, selected_species] .= 1
    end

    return psquads

end