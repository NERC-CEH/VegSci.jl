# Fidelity Measures


``` julia
using VegSci
using Clustering
using CSV
using DataFrames
using NamedArrays
```

    [ Info: Precompiling VegSci [31ef3d96-7648-423b-b844-5d3d759e611e]
    [ Info: Precompiling Clustering [aaaa29a8-35af-508c-8bc3-b662a17a0fe5]

## Outline

This methods page provides an overview of the step-by-step
implementation of the fidelity metrics following Chytrý et al. (2002).

The functionality outlined here aims to replicate the functionality
found in the program JUICE Tichý (2002).

## Load example data.

Load the `dune` dataset bundled with the R package `vegan`.

``` julia
begin
    dune_df = CSV.read("../data/dune.csv", DataFrame, header = 1)
    dune_na = NamedArray(Array(dune_df))
    NamedArrays.setnames!(dune_na, names(dune_df), 2)
  dune_na = float.(dune_na)
end
```

    20×30 Named Matrix{Float64}
    A ╲ B │ Achimill  Agrostol  Airaprae  …  Vicilath  Bracruta  Callcusp
    ──────┼──────────────────────────────────────────────────────────────
    1     │      1.0       0.0       0.0  …       0.0       0.0       0.0
    2     │      3.0       0.0       0.0          0.0       0.0       0.0
    3     │      0.0       4.0       0.0          0.0       2.0       0.0
    4     │      0.0       8.0       0.0          0.0       2.0       0.0
    5     │      2.0       0.0       0.0          0.0       2.0       0.0
    6     │      2.0       0.0       0.0          0.0       6.0       0.0
    7     │      2.0       0.0       0.0          0.0       2.0       0.0
    8     │      0.0       4.0       0.0          0.0       2.0       0.0
    9     │      0.0       3.0       0.0          0.0       2.0       0.0
    10    │      4.0       0.0       0.0          1.0       2.0       0.0
    11    │      0.0       0.0       0.0          2.0       4.0       0.0
    12    │      0.0       4.0       0.0          0.0       4.0       0.0
    13    │      0.0       5.0       0.0          0.0       0.0       0.0
    14    │      0.0       4.0       0.0          0.0       0.0       4.0
    15    │      0.0       4.0       0.0          0.0       4.0       0.0
    16    │      0.0       7.0       0.0          0.0       4.0       3.0
    17    │      2.0       0.0       2.0          0.0       0.0       0.0
    18    │      0.0       0.0       0.0          1.0       6.0       0.0
    19    │      0.0       0.0       3.0          0.0       3.0       0.0
    20    │      0.0       5.0       0.0  …       0.0       4.0       3.0

## Identify Clusters

Let’s identify some clusters, storing the cluster-releve memberships as
a dictionary.

``` julia
r = Clustering.fuzzy_cmeans(transpose(dune_na), 3, 2)

cluster_weights = r.weights
clusters_vec = vec(Tuple.(findmax(cluster_weights, dims = 2)[2]))
clusters_mat = hcat(first.(clusters_vec), last.(clusters_vec))

clusters = Dict

for i in unique(clusters_mat[:,2])

    rowids = clusters_mat[clusters_mat[:,2] .== i, :][:,1]
    rownames = names(dune_na)[1][rowids]
    clusters_i = Dict(string(i) => string.(rownames))
    clusters = merge(clusters, clusters_i)
    
end

clusters
```

    Dict{String, Vector{String}} with 3 entries:
      "1" => ["3", "4", "8", "9", "12", "13"]
      "2" => ["14", "15", "16", "19", "20"]
      "3" => ["1", "2", "5", "6", "7", "10", "11", "17", "18"]

Create a presence-absence matrix by replacing all non-zero values with
an Integer value of 1.

``` julia
dune_pa = Int.(dune_na .!= 0)
```

    20×30 Named Matrix{Int64}
    A ╲ B │ Achimill  Agrostol  Airaprae  …  Vicilath  Bracruta  Callcusp
    ──────┼──────────────────────────────────────────────────────────────
    1     │        1         0         0  …         0         0         0
    2     │        1         0         0            0         0         0
    3     │        0         1         0            0         1         0
    4     │        0         1         0            0         1         0
    5     │        1         0         0            0         1         0
    6     │        1         0         0            0         1         0
    7     │        1         0         0            0         1         0
    8     │        0         1         0            0         1         0
    9     │        0         1         0            0         1         0
    10    │        1         0         0            1         1         0
    11    │        0         0         0            1         1         0
    12    │        0         1         0            0         1         0
    13    │        0         1         0            0         0         0
    14    │        0         1         0            0         0         1
    15    │        0         1         0            0         1         0
    16    │        0         1         0            0         1         1
    17    │        1         0         1            0         0         0
    18    │        0         0         0            1         1         0
    19    │        0         0         1            0         1         0
    20    │        0         1         0  …         0         1         1

## Calculate fidelity metrics

Using cluster 1 as an example, *p* = 1.

### Observed Frequencies

First, calculate the observed frequencies *f*(*o*)<sub>*i*</sub>.

***N***

Calculate *N* the number of releves in the matrix.

``` julia
N = size(dune_pa)[1]
```

    20

***N*<sub>*p*</sub>**

Calculate *N*<sub>*p*</sub> the number of releves in the particular
vegetation unit.

``` julia
Np = length(getindex(clusters, "1"))
```

    6

***n***

Calculate *n* the number of occurrences of the species in the matrix.

``` julia
n = sum(dune_pa, dims = 1)
setnames!(n, ["all"], 1)
n
```

    1×30 Named Matrix{Int64}
    A ╲ B │ Achimill  Agrostol  Airaprae  …  Vicilath  Bracruta  Callcusp
    ──────┼──────────────────────────────────────────────────────────────
    all   │        7        10         2  …         3        15         3

***f*(*o*)<sub>1</sub>**

Calculate *f*(*o*)<sub>1</sub> or *n*<sub>*p*</sub>, the number of
occurences of the species in the particular vegetation unit.

``` julia
np = sum(dune_pa[getindex(clusters, "1"),:], dims = 1)
fo_1 = np
```

    1×30 Named Matrix{Int64}
     A ╲ B │ Achimill  Agrostol  Airaprae  …  Vicilath  Bracruta  Callcusp
    ───────┼──────────────────────────────────────────────────────────────
    sum(A) │        0         6         0  …         0         5         0

### Dufrêne-Legendre Indicator Value Index

Calculate the Dufrêne-Legendre Indicator Value Index (IndVal) measure of
species fidelity (Dufrêne and Legendre 1997) which does not use observed
and expected frequencies.

``` julia
indval = ((np .* (N - Np)) ./ (((n .* Np) .- (2 .* np)) .+ (np .* N))) .* (np ./ Np)
```

    1×30 Named Matrix{Float64}
     A ╲ B │  Achimill   Agrostol   Airaprae  …   Vicilath   Bracruta   Callcusp
    ───────┼────────────────────────────────────────────────────────────────────
    sum(A) │       0.0        0.5        0.0  …        0.0   0.324074        0.0

***f*(*o*)<sub>2</sub>**

Calculate *n* − *n*<sub>*p*</sub> the number of releves containing the
species that aren’t in a particular vegetation unit.

``` julia
fo_2 = n - np
```

    ┌ Warning: Using names of left argument
    └ @ NamedArrays ~/.julia/packages/NamedArrays/nqYoX/src/arithmetic.jl:25

    1×30 Named Matrix{Int64}
    A ╲ B │ Achimill  Agrostol  Airaprae  …  Vicilath  Bracruta  Callcusp
    ──────┼──────────────────────────────────────────────────────────────
    all   │        7         4         2  …         3        10         3

***f*(*o*)<sub>3</sub>**

Calculate *N*<sub>*p*</sub> − *n*<sub>*p*</sub> the number of releves
which don’t contain the species in a particular vegetation unit.

``` julia
fo_3 = Np .- np
```

    1×30 Named Matrix{Int64}
     A ╲ B │ Achimill  Agrostol  Airaprae  …  Vicilath  Bracruta  Callcusp
    ───────┼──────────────────────────────────────────────────────────────
    sum(A) │        6         0         6  …         6         1         6

***f*(*o*)<sub>4</sub>**

Calculate *N* − *N*<sub>*p*</sub> − *n* + *n*<sub>*p*</sub> the number
of releves not containing the species and not in the vegetation unit

``` julia
fo_4 = N - Np .- n .+ np
```

    1×30 Named Matrix{Int64}
    A ╲ B │ Achimill  Agrostol  Airaprae  …  Vicilath  Bracruta  Callcusp
    ──────┼──────────────────────────────────────────────────────────────
    all   │        7        10        12  …        11         4        11

### Expected Frequencies

Then calculate the expected frequencies *f*(*e*)<sub>*i*</sub>

\*\*$f(e)\_{1}\*\*

Calculate the expected number of the releves containing the species in
the particular vegetation unit *n* ⋅ *N*<sub>*p*</sub>/*N*.

``` julia
fe_1 = n .* (Np / N)
```

    1×30 Named Matrix{Float64}
    A ╲ B │ Achimill  Agrostol  Airaprae  …  Vicilath  Bracruta  Callcusp
    ──────┼──────────────────────────────────────────────────────────────
    all   │      2.1       3.0       0.6  …       0.9       4.5       0.9

\*\*$f(e)\_{2}\*\*

Cal

``` julia
fe_2 = n .* ((N - Np) / N)
```

    1×30 Named Matrix{Float64}
    A ╲ B │ Achimill  Agrostol  Airaprae  …  Vicilath  Bracruta  Callcusp
    ──────┼──────────────────────────────────────────────────────────────
    all   │      4.9       7.0       1.4  …       2.1      10.5       2.1

\*\*$f(e)\_{3}\*\*

``` julia
fe_3 = (N .- n) .* (Np ./ N)
```

    1×30 Named Matrix{Float64}
    A ╲ B │ Achimill  Agrostol  Airaprae  …  Vicilath  Bracruta  Callcusp
    ──────┼──────────────────────────────────────────────────────────────
    all   │      3.9       3.0       5.4  …       5.1       1.5       5.1

\*\*$f(e)\_{4}\*\*

``` julia
fe_4 = (N .- n) .* ((N - Np) / N)
```

    1×30 Named Matrix{Float64}
    A ╲ B │ Achimill  Agrostol  Airaprae  …  Vicilath  Bracruta  Callcusp
    ──────┼──────────────────────────────────────────────────────────────
    all   │      9.1       7.0      12.6  …      11.9       3.5      11.9

## Compose functions

Here we create two functions, the first which doesn’t use observed and
expected frequencies - the Dufrêne-Legendre Indicator Value Index
function.

``` julia
function indval_fidelity(x::NamedMatrix, clusters::Dict{String, Vector{String}})

    N = size(x)[1]
    n = sum(x, dims = 1)

    indval_all = NamedArrays.NamedArray(zeros(length(names(clusters)), size(x)[2]), names = (names(clusters), names(x)[2]))

    for i in names(clusters)
        Np = length(getindex(clusters, i))
        np = sum(x[getindex(clusters, i),:], dims = 1)
        indval = ((np .* (N - Np)) ./ (((n .* Np) .- (2 .* np)) .+ (np .* N))) .* (np ./ Np)
        indval_all[i,:] = indval
    end

    return indval_all

end
```

    indval_fidelity (generic function with 1 method)

Test the `indval_fidelity` function.

``` julia
indval_fidelity(dune_pa, clusters)
```

    3×30 Named Matrix{Float64}
    A ╲ B │  Achimill   Agrostol   Airaprae  …   Vicilath   Bracruta   Callcusp
    ──────┼────────────────────────────────────────────────────────────────────
    1     │       0.0        0.5        0.0  …        0.0   0.324074        0.0
    2     │       0.0   0.393443   0.107143           0.0   0.326531   0.391304
    3     │  0.316872        0.0  0.0339506  …   0.135802    0.18107        0.0

## References

Chytrý, Milan, Lubomír Tichý, Jason Holt, and Zoltán Botta-Dukát. 2002.
“Determination of Diagnostic Species with Statistical Fidelity
Measures.” *Journal of Vegetation Science* 13 (1): 79–90.
<https://doi.org/10.1111/j.1654-1103.2002.tb02025.x>.

Dufrêne, Marc, and Pierre Legendre. 1997. “Species Assemblages and
Indicator Species:the Need for a Flexible Asymmetrical Approach.”
*Ecological Monographs* 67 (3): 345–66.
[https://doi.org/10.1890/0012-9615(1997)067\[0345:SAAIST\]2.0.CO;2](https://doi.org/10.1890/0012-9615(1997)067[0345:SAAIST]2.0.CO;2).

Tichý, Lubomír. 2002. “JUICE, Software for Vegetation Classification.”
*Journal of Vegetation Science* 13 (3): 451–53.
<https://doi.org/10.1111/j.1654-1103.2002.tb02069.x>.
