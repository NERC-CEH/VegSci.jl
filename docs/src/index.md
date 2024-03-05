

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
x = VegSci.generate_test_array(rown = 30, coln = 20, meancoloccs = 10, rowprefix = "SiteA-", colprefix = "Species")
```

    30×20 Named Matrix{Float64}
    Releve ╲ Species │    Species1     Species2  …    Species19    Species20
    ─────────────────┼──────────────────────────────────────────────────────
    SiteA-1          │         0.0     0.130287  …     0.104092          0.0
    SiteA-2          │         0.0    0.0575279        0.126437          0.0
    SiteA-3          │    0.135692     0.126111             0.0          0.0
    SiteA-4          │         0.0     0.228196             0.0          0.0
    SiteA-5          │         0.0          0.0             0.0     0.129159
    SiteA-6          │         0.0    0.0696195        0.159017    0.0175125
    SiteA-7          │         0.0    0.0100315             0.0     0.107451
    SiteA-8          │    0.123417    0.0887687             0.0          0.0
    SiteA-9          │         0.0     0.123951        0.110472    0.0810747
    SiteA-10         │    0.170721    0.0522378        0.148841    0.0838056
    SiteA-11         │         0.0    0.0381516             0.0     0.151037
    ⋮                            ⋮            ⋮  ⋱            ⋮            ⋮
    SiteA-20         │         0.0    0.0645953             0.0     0.145417
    SiteA-21         │    0.157137          0.0             0.0    0.0812415
    SiteA-22         │   0.0251965          0.0        0.152538          0.0
    SiteA-23         │         0.0     0.284651             0.0          0.0
    SiteA-24         │   0.0192188     0.115009             0.0    0.0657227
    SiteA-25         │         0.0     0.137841       0.0377793     0.153853
    SiteA-26         │         0.0     0.107912             0.0          0.0
    SiteA-27         │         0.0     0.139419             0.0          0.0
    SiteA-28         │   0.0373618       0.0713             0.0    0.0592767
    SiteA-29         │    0.156715    0.0565809             0.0    0.0503391
    SiteA-30         │   0.0664498          0.0  …     0.174712          0.0

### Classification

Let’s identify some clusters.

``` julia
r = Clustering.fuzzy_cmeans(transpose(x), 3, 2)

cluster_weights = r.weights
memberships_vec = vec(Tuple.(findmax(cluster_weights, dims = 2)[2]))
memberships_mat = hcat(first.(memberships_vec), last.(memberships_vec))

memberships = Dict

for i in unique(memberships_mat[:,2])

    rowids = memberships_mat[memberships_mat[:,2] .== i, :][:,1]
    memberships_i = Dict(i => rowids)
    memberships = merge(memberships, memberships_i)
    
end

memberships
```

    Dict{Int64, Vector{Int64}} with 3 entries:
      2 => [1, 4, 7, 9, 11, 18, 20, 23, 25, 26, 28, 30]
      3 => [2, 3, 6, 8, 13, 16, 17, 21, 24, 27, 29]
      1 => [5, 10, 12, 14, 15, 19, 22]

### Creation of Syntopic Tables

Once the plots have been grouped into clusters, we can proceed to
summarise their composition via the creation of `SyntopicTable` objects.

``` julia
syn_1 = VegSci.compose_syntopic_table_object("Syn1", x[getindex(memberships, 1),:])
syn_2 = VegSci.compose_syntopic_table_object("Syn2", x[getindex(memberships, 2),:])
VegSci.print_summary_syntopic_table(syn_2, "normal", "cover_proportion")
```



    Community Name: Syn2
    Releves: n = 12
    Species: n = 20
    ┌───────────┬───────────────────┬─────────────────┐
    │   Species │ RelativeFrequency │       Abundance │
    ├───────────┼───────────────────┼─────────────────┤
    │  Species2 │               0.8 │ 0.1 (0.0 - 0.3) │
    │ Species12 │               0.8 │ 0.1 (0.0 - 0.2) │
    │ Species14 │               0.8 │ 0.1 (0.1 - 0.2) │
    │ Species17 │               0.8 │ 0.1 (0.0 - 0.2) │
    │  Species9 │               0.7 │ 0.1 (0.0 - 0.1) │
    │ Species10 │               0.7 │ 0.1 (0.0 - 0.2) │
    │ Species13 │               0.7 │ 0.1 (0.0 - 0.2) │
    │ Species18 │               0.7 │ 0.1 (0.0 - 0.3) │
    │  Species4 │               0.6 │ 0.1 (0.0 - 0.2) │
    │  Species7 │               0.5 │ 0.1 (0.0 - 0.2) │
    │  Species8 │               0.5 │ 0.1 (0.0 - 0.2) │
    │ Species20 │               0.5 │ 0.1 (0.1 - 0.2) │
    │  Species3 │               0.4 │ 0.1 (0.0 - 0.2) │
    │  Species6 │               0.3 │ 0.1 (0.1 - 0.1) │
    │ Species15 │               0.3 │ 0.1 (0.0 - 0.3) │
    │ Species19 │               0.3 │ 0.1 (0.0 - 0.2) │
    │  Species5 │               0.3 │ 0.0 (0.0 - 0.0) │
    │ Species11 │               0.2 │ 0.1 (0.0 - 0.1) │
    │  Species1 │               0.2 │ 0.0 (0.0 - 0.1) │
    │ Species16 │               0.1 │ 0.1 (0.1 - 0.1) │
    └───────────┴───────────────────┴─────────────────┘

### Identification of High-Fidelity Species

### Generation of Pseudo-Releves

### Assignment of Releves to Vegetation Classes

Let’s generate a second example matrix, consisting of sample 5 releves,
against which we want to calculate the similarity.

``` julia
y = VegSci.generate_test_array(rown = 5, coln = 30, meancoloccs = 5, rowprefix = "SiteB-", colprefix = "Species")
```

    5×30 Named Matrix{Float64}
    Releve ╲ Species │   Species1    Species2  …   Species29   Species30
    ─────────────────┼──────────────────────────────────────────────────
    SiteB-1          │  0.0297561         0.0  …         0.0         0.0
    SiteB-2          │        0.0    0.173407            0.0         0.0
    SiteB-3          │        0.0         0.0      0.0633288         0.0
    SiteB-4          │        0.0         0.0            0.0         0.0
    SiteB-5          │        0.0         0.0  …   0.0526365    0.151109

Three methods will be demonstrated.

### Jaccard Similarity

### Steinhaus coefficient

First, let’s compose a syntopic table object from the “y” sample data
and extract the syntopic tables in matrix format.

``` julia
syn_y = VegSci.compose_syntopic_table_object("Sample", y)
syn_y_mat = VegSci.extract_syntopic_matrix(syn_y)
syn_1_mat = VegSci.extract_syntopic_matrix(syn_1)
syn_2_mat = VegSci.extract_syntopic_matrix(syn_2)
```

    1×20 Named Matrix{Float64}
    A ╲ B │  Species1   Species2   Species3  …  Species18  Species19  Species20
    ──────┼────────────────────────────────────────────────────────────────────
    Syn2  │ 0.0472968   0.115932  0.0741902  …    0.10908   0.107282   0.126434

Now we have three matrices, containg the relative frequencies of each
species present in the sample releves which constitute each syntaxon.
However, each of the syntaxa are composed of a different set of species,
in Julia we need a helper function to merge these matrices and ensure
each matrix contains each species across all the matrices. This function
is broadly equivalent to the R function `base::merge`.

``` julia
merged_syn_mats = VegSci.merge_namedarrays([syn_y_mat, syn_1_mat, syn_2_mat])
```

    3×27 Named Matrix{Float64}
     A ╲ B │   Species1    Species2  …   Species19   Species12
    ───────┼──────────────────────────────────────────────────
    Sample │  0.0297561    0.173407  …         0.0         0.0
    Syn1   │   0.145251   0.0794433       0.127773         0.0
    Syn2   │  0.0472968    0.115932  …    0.107282   0.0904342

``` julia
VegSci.steinhaus_coefficient(merged_syn_mats[[:"Sample"],:], merged_syn_mats[Not(:"Sample"), :])
```

    1×2 Named Matrix{Float64}
     A ╲ B │     Syn1      Syn2
    ───────┼───────────────────
    Sample │ 0.337378  0.324535

### Multivariate Analysis

### Ecological Trajectory Analysis

### Example Workflow

## External Resources

## Implemented Methodologies

## Contribute

## Acknowledgements

## References

Bezanson, Jeff, Alan Edelman, Stefan Karpinksi, and Viral B. Shah. 2017.
“Julia: A Fresh Approach to Numerical Computing.” *SIAM Review* 59 (1):
65–98. <https://doi.org/10.1137/141000671>.

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
