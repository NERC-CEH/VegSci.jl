"""
    indval_fidelity(x::NamedMatrix, clusters::Dict{String, Vector{String}})::NamedMatrix

Calculate the Dufrêne-Legendre Indicator Value Index (IndVal) fidelity values for a given releve by species matrix
and set of clusters.

### Input

- `x` -- A releve by species matrix of the class NamedArrays::NamedMatrix
- `clusters` -- A dictionary containing the names of the clusters and associated lists of row names belonging to that cluster of the class Dict{String, Vector{String}}

### Output

A cluster by species NamedMatrix containing the fidelity values.

### Notes


### Algorithm

This function implements the Dufrêne-Legendre Indicator Value Index as described in Dufrêne & Legendre (1997)

### References

Dufrêne, M. & Legendre, P. 1997. Species assemblages and indicator species: the need for a flexible asymmetrical approach. Ecol. Monogr. 67: 345-366.

"""
function indval_fidelity(x::NamedMatrix, clusters::Dict{String, Vector{String}})

    N = size(x)[1]
    n = sum(x, dims = 1)

    indval_all = NamedArrays.NamedArray(zeros(length(names(clusters)), size(x)[2]), names = (names(clusters), names(x)[2]))

    for i in names(clusters)
        Np = length(getindex(clusters, i))
        np = sum(x[getindex(clusters, i),:], dims = 1)
        indval = ((np .* (N .- Np)) ./ (((n .* Np) .- (2 .* np)) .+ (np .* N))) .* (np ./ Np)
        indval_all[i,:] = indval
    end

    return indval_all

end