

``` @meta
CurrentModule = VegSci
```

Documentation for
[VegSci.jl](https://github.com/ZekeMarhshall/VegSci.jl).

# VegSci

Tools for vegetation science.

## Background

`VegSci.jl` is a package containing tools for vegetation science using
Julia (Bezanson et al. 2017), a growing scientific programming language
which solves the ‘two language problem’ (Roesch et al. 2023), offering C
and FORTRAN-like performance alongside the readability and
user-friendliness of higher level languages such as Python. `VegSci.jl`
aims to collate selected functionality found in popular vegetation
science software programs/packages such as JUICE, vegan, ade4, vegclust,
vegsoup, and ecotraj into a single location with a user-friendly API and
transparent methodologies. `VegSci.jl` is being developed with the aim
of assisting in the creation of high-performance, reproducible
analytical pipelines in vegetation research (Sperandii et al. 2024),
developed primarily with the application to the vegetation of Great
Britain in mind, but fully generalisable. Nomenclature follows
Theurillat et al. (2021).

## Installation

To install the latest stable release of `VegSci`:

``` julia
using Pkg
Pkg.add("VegSci")
```

To install the development version of `VegSci`:

``` julia
using Pkg
Pkg.add(url="https://github.com/ZekeMarshall/VegSci.jl", rev = "develop")
```

## Usage Example

To demonstrate…

First we begin with generating two example plot by species
`NamedArrays.NamedMatrix` object using the function
`VegSci.generate_test_array` as test data.

``` julia
dune = CSV.read("./data/dune.csv", DataFrame, header = 1)
x = float(NamedArray(Array(dune)))
NamedArrays.setnames!(x, names(dune), 2)
x
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

### Classification

Let’s identify some potential associations using fuzzy c-means
clustering, storing the cluster-releve memberships as a dictionary.

``` julia
r = Clustering.fuzzy_cmeans(transpose(x), 3, 2)

cluster_weights = r.weights
clusters_vec = vec(Tuple.(findmax(cluster_weights, dims = 2)[2]))
clusters_mat = hcat(first.(clusters_vec), last.(clusters_vec))

clusters = Dict

for i in unique(clusters_mat[:,2])

    rowids = clusters_mat[clusters_mat[:,2] .== i, :][:,1]
    rownames = names(x)[1][rowids]
    clusters_i = Dict(string(i) => string.(rownames))
    clusters = merge(clusters, clusters_i)
    
end

clusters
```

    Dict{String, Vector{String}} with 3 entries:
      "1" => ["3", "4", "8", "9", "12", "13"]
      "2" => ["1", "2", "5", "6", "7", "10", "11", "17", "18"]
      "3" => ["14", "15", "16", "19", "20"]

### Creation of Syntopic Tables

Once the plots have been grouped into clusters, we can proceed to
summarise their composition via the creation of `SyntopicTable` objects.

``` julia
syn_1 = VegSci.compose_syntopic_table_object("Syn1", x[getindex(clusters, "1"),:])
syn_2 = VegSci.compose_syntopic_table_object("Syn2", x[getindex(clusters, "2"),:])
VegSci.print_summary_syntopic_table(syn_1, "normal", "proportion")
```



    Community Name: Syn1
    Releves: n = 6
    Species: n = 19
    Postive Indicators: 
    Negative Indicators: 
    ┌──────────┬───────────────────┬─────────────────┐
    │  Species │ RelativeFrequency │       Abundance │
    ├──────────┼───────────────────┼─────────────────┤
    │ Alopgeni │               1.0 │ 5.0 (2.0 - 8.0) │
    │  Poatriv │               1.0 │ 5.0 (4.0 - 9.0) │
    │ Agrostol │               1.0 │ 4.0 (3.0 - 8.0) │
    │ Scorautu │               1.0 │ 2.0 (2.0 - 3.0) │
    │ Trifrepe │               1.0 │ 2.0 (1.0 - 3.0) │
    │  Poaprat │               0.8 │ 4.0 (2.0 - 5.0) │
    │ Sagiproc │               0.8 │ 2.0 (2.0 - 5.0) │
    │ Bracruta │               0.8 │ 2.0 (2.0 - 4.0) │
    │ Lolipere │               0.7 │ 4.5 (2.0 - 6.0) │
    │ Elymrepe │               0.5 │ 4.0 (4.0 - 6.0) │
    │ Juncbufo │               0.5 │ 4.0 (3.0 - 4.0) │
    │ Juncarti │               0.3 │ 4.0 (4.0 - 4.0) │
    │ Bellpere │               0.3 │ 2.0 (2.0 - 2.0) │
    │ Ranuflam │               0.3 │ 2.0 (2.0 - 2.0) │
    │ Rumeacet │               0.3 │ 2.0 (2.0 - 2.0) │
    │ Eleopalu │               0.2 │ 4.0 (4.0 - 4.0) │
    │ Bromhord │               0.2 │ 3.0 (3.0 - 3.0) │
    │ Cirsarve │               0.2 │ 2.0 (2.0 - 2.0) │
    │ Chenalbu │               0.2 │ 1.0 (1.0 - 1.0) │
    └──────────┴───────────────────┴─────────────────┘

### Species Fidelity

Calculate the species fidelity scores using the Dufrêne-Legendre
Indicator Value Index (IndVal) metric for each cluster.

``` julia
indval_fidelity = VegSci.indval_fidelity(x, clusters)
```

    3×30 Named Matrix{Float64}
    A ╲ B │  Achimill   Agrostol   Airaprae  …   Vicilath   Bracruta   Callcusp
    ──────┼────────────────────────────────────────────────────────────────────
    1     │       0.0    2.30976        0.0  …        0.0   0.658824        0.0
    2     │   0.72428        0.0  0.0603567       0.18107   0.706757        0.0
    3     │       0.0        2.0   0.341772  …        0.0    1.31068    1.30435

Let’s update the syntopic table objects with the fidelity values. By
default the IndVal fidelity measure does not distinguish negative
fidelity (Chytrý et al. 2002), so we just provide a ‘cut value’ to
extract the postive, high-fidelity species.

``` julia
syn_1_f_mat = indval_fidelity[["1"], :]
p_cut = 0.5
syn_1.fidelity = vec(syn_1_f_mat)
syn_1.fidelity_p = names(syn_1_f_mat[:, vec(map(col -> any(col .>= 0.5), eachcol(syn_1_f_mat)))])[2]
VegSci.print_summary_syntopic_table(syn_1, "normal", "proportion")
```



    Community Name: Syn1
    Releves: n = 6
    Species: n = 19
    Postive Indicators: Agrostol Alopgeni Elymrepe Juncarti Juncbufo Lolipere Poaprat Poatriv Sagiproc Scorautu Trifrepe Bracruta
    Negative Indicators: 
    ┌──────────┬───────────────────┬─────────────────┐
    │  Species │ RelativeFrequency │       Abundance │
    ├──────────┼───────────────────┼─────────────────┤
    │ Alopgeni │               1.0 │ 5.0 (2.0 - 8.0) │
    │  Poatriv │               1.0 │ 5.0 (4.0 - 9.0) │
    │ Agrostol │               1.0 │ 4.0 (3.0 - 8.0) │
    │ Scorautu │               1.0 │ 2.0 (2.0 - 3.0) │
    │ Trifrepe │               1.0 │ 2.0 (1.0 - 3.0) │
    │  Poaprat │               0.8 │ 4.0 (2.0 - 5.0) │
    │ Sagiproc │               0.8 │ 2.0 (2.0 - 5.0) │
    │ Bracruta │               0.8 │ 2.0 (2.0 - 4.0) │
    │ Lolipere │               0.7 │ 4.5 (2.0 - 6.0) │
    │ Elymrepe │               0.5 │ 4.0 (4.0 - 6.0) │
    │ Juncbufo │               0.5 │ 4.0 (3.0 - 4.0) │
    │ Juncarti │               0.3 │ 4.0 (4.0 - 4.0) │
    │ Bellpere │               0.3 │ 2.0 (2.0 - 2.0) │
    │ Ranuflam │               0.3 │ 2.0 (2.0 - 2.0) │
    │ Rumeacet │               0.3 │ 2.0 (2.0 - 2.0) │
    │ Eleopalu │               0.2 │ 4.0 (4.0 - 4.0) │
    │ Bromhord │               0.2 │ 3.0 (3.0 - 3.0) │
    │ Cirsarve │               0.2 │ 2.0 (2.0 - 2.0) │
    │ Chenalbu │               0.2 │ 1.0 (1.0 - 1.0) │
    └──────────┴───────────────────┴─────────────────┘

### Generation of Pseudo-Releves

### Assignment of Releves to Vegetation Classes

Let’s compare the similarity of `syn_1` and `syn_2`.

### Steinhaus coefficient

First, let’s compose a syntopic table object from the “y” sample data
and extract the syntopic tables in matrix format.

``` julia
syn_1_mat = VegSci.extract_syntopic_matrix(syn_1)
syn_2_mat = VegSci.extract_syntopic_matrix(syn_2)
```

    1×21 Named Matrix{Float64}
    A ╲ B │ Achimill  Airaprae  Alopgeni  …  Trifrepe  Vicilath  Bracruta
    ──────┼──────────────────────────────────────────────────────────────
    Syn2  │      2.0       2.0       2.0  …       3.0       1.0       3.0

Now we have three matrices, containing the relative frequencies of each
species present in the sample releves which constitute each syntaxon.
However, each of the syntaxa are composed of a different set of species,
in Julia we need a helper function to merge these matrices and ensure
each matrix contains each species across all the matrices. This function
is broadly equivalent to the R function `base::merge`.

``` julia
merged_syn_mats = VegSci.merge_namedarrays([syn_1_mat, syn_2_mat])
```

    2×27 Named Matrix{Float64}
    A ╲ B │ Agrostol  Alopgeni  Bellpere  …  Salirepe  Trifprat  Vicilath
    ──────┼──────────────────────────────────────────────────────────────
    Syn1  │      4.0       5.0       2.0  …       0.0       0.0       0.0
    Syn2  │      0.0       2.0       2.0  …       3.0       2.0       1.0

``` julia
VegSci.steinhaus_coefficient(merged_syn_mats[[:"Syn1"], :], merged_syn_mats[Not(:"Syn1"), :])
```

    1×1 Named Matrix{Float64}
    A ╲ B │     Syn2
    ──────┼─────────
    Syn1  │ 0.595041

### Multivariate Analysis

### Example Workflow

## External Resources

## Implemented Methodologies

## Contribute

## Acknowledgements

## References

Bezanson, Jeff, Alan Edelman, Stefan Karpinksi, and Viral B. Shah. 2017.
“Julia: A Fresh Approach to Numerical Computing.” *SIAM Review* 59 (1):
65–98. <https://doi.org/10.1137/141000671>.

Chytrý, Milan, Lubomír Tichý, Jason Holt, and Zoltán Botta-Dukát. 2002.
“Determination of Diagnostic Species with Statistical Fidelity
Measures.” *Journal of Vegetation Science* 13 (1): 79–90.
<https://doi.org/10.1111/j.1654-1103.2002.tb02025.x>.

Roesch, Elisabeth, Joe G. Greener, Adam L. MacLean, Huda Nassar,
Christopher Rackauckas, Timothy E. Holy, and Michael P. H. Stumpf. 2023.
“Julia for Biologists.” *Nature Methods* 20 (5): 655–64.
<https://doi.org/10.1038/s41592-023-01832-z>.

Sperandii, Marta Gaia, Manuele Bazzichetto, Glenda Mendieta-Leiva,
Sebastian Schmidtlein, Michael Bott, Renato A. Ferreira de Lima, Valério
D. Pillar, Jodi N. Price, Viktoria Wagner, and Milan Chytrý. 2024.
“Towards More Reproducibility in Vegetation Research.” *Journal of
Vegetation Science* 35 (1): e13224. <https://doi.org/10.1111/jvs.13224>.

Theurillat, Jean-Paul, Wolfgang Willner, Federico Fernández-González,
Helga Bültmann, Andraž Čarni, Daniela Gigante, Ladislav Mucina, and
Heinrich Weber. 2021. “International Code of Phytosociological
Nomenclature. 4th Edition.” *Applied Vegetation Science* 24 (1): e12491.
<https://doi.org/10.1111/avsc.12491>.
