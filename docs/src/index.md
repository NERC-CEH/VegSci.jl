

``` @meta
CurrentModule = VegSci
```

Documentation for
[VegSci.jl](https://github.com/ZekeMarhshall/VegSci.jl).

# VegSci

Tools for vegetation science.

## Background

`VegSci.jl` is a package containing tools for vegetation science using
the Julia (Bezanson et al. 2017), a growing scientific programming
language which solves the ‘two language problem’ (@ Roesch et al. 2023),
offering C and FORTRAN-like performance alongside the readability and
user-friendliness of higher level languages such as Python. `VegSci.jl`
aims to collate selected functionality found in popular vegetation
science software programs/packages such as JUICE, vegan, ade4, vegclust,
vegsoup, and ecotraj into a single location with a user-friendly API and
transparent methodologies. `VegSci.jl` is being developed with the aim
of assisting in the creation of high-performance, reproducible analysis
pipelines in vegetation research (Sperandii et al. 2024), developed
primarily with the application to the vegetation of Great Britain in
mind, but fully generalisable. Nomenclature follows Theurillat et al.
(2021).

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
    Releve ╲ Species │   Species1    Species2  …   Species19   Species20
    ─────────────────┼──────────────────────────────────────────────────
    SiteA-1          │        0.0   0.0961665  …         0.0         0.0
    SiteA-2          │        0.0         0.0       0.101352         0.0
    SiteA-3          │   0.232736         0.0            0.0         0.0
    SiteA-4          │        0.0         0.0      0.0112884    0.113199
    SiteA-5          │  0.0343561    0.151313       0.119557         0.0
    SiteA-6          │  0.0482408         0.0       0.132966   0.0459301
    SiteA-7          │  0.0481567         0.0      0.0787806   0.0791549
    SiteA-8          │        0.0         0.0       0.163913         0.0
    SiteA-9          │   0.119563   0.0542965            0.0   0.0582056
    SiteA-10         │        0.0  0.00102027            0.0   0.0605912
    SiteA-11         │        0.0         0.0       0.194125         0.0
    ⋮                           ⋮           ⋮  ⋱           ⋮           ⋮
    SiteA-20         │   0.157566   0.0929004      0.0696219  0.00322057
    SiteA-21         │   0.149004   0.0849708            0.0   0.0396169
    SiteA-22         │        0.0   0.0806118      0.0954818    0.199833
    SiteA-23         │        0.0         0.0       0.248834   0.0484612
    SiteA-24         │   0.164072         0.0      0.0217969         0.0
    SiteA-25         │  0.0938849   0.0388366            0.0   0.0323463
    SiteA-26         │     0.1296   0.0720648            0.0         0.0
    SiteA-27         │   0.111855    0.153854       0.153939   0.0567083
    SiteA-28         │   0.100784         0.0            0.0         0.0
    SiteA-29         │        0.0         0.0            0.0         0.0
    SiteA-30         │    0.14317   0.0628497  …         0.0    0.144168

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
      2 => [18, 19, 27]
      3 => [3, 9, 10, 12, 13, 14, 16, 20, 21, 24, 25, 26, 28, 30]
      1 => [1, 2, 4, 5, 6, 7, 8, 11, 15, 17, 22, 23, 29]

### Creation of Syntopic Tables

Once the plots have been grouped into clusters, we can proceed to
summarise their composition via the creation of `SyntopicTable` objects.

``` julia
syn_1 = VegSci.compose_syntopic_table_object("Syn1", x[getindex(memberships, 1),:])
syn_2 = VegSci.compose_syntopic_table_object("Syn2", x[getindex(memberships, 2),:])
VegSci.print_summary_syntopic_table(syn_2, "normal", "cover_proportion")
```



    Community Name: Syn2
    Releves: n = 3
    Species: n = 17
    ┌───────────┬───────────────────┬─────────────────┐
    │   Species │ RelativeFrequency │       Abundance │
    ├───────────┼───────────────────┼─────────────────┤
    │  Species1 │               1.0 │ 0.1 (0.0 - 0.2) │
    │ Species18 │               1.0 │ 0.1 (0.1 - 0.2) │
    │  Species7 │               0.7 │ 0.2 (0.2 - 0.2) │
    │  Species5 │               0.7 │ 0.1 (0.1 - 0.1) │
    │  Species8 │               0.7 │ 0.1 (0.0 - 0.2) │
    │ Species10 │               0.7 │ 0.1 (0.0 - 0.1) │
    │ Species13 │               0.7 │ 0.1 (0.0 - 0.2) │
    │ Species15 │               0.7 │ 0.1 (0.0 - 0.2) │
    │ Species17 │               0.7 │ 0.1 (0.1 - 0.1) │
    │ Species19 │               0.7 │ 0.1 (0.0 - 0.2) │
    │  Species2 │               0.3 │ 0.2 (0.2 - 0.2) │
    │  Species6 │               0.3 │ 0.2 (0.2 - 0.2) │
    │  Species3 │               0.3 │ 0.1 (0.1 - 0.1) │
    │ Species20 │               0.3 │ 0.1 (0.1 - 0.1) │
    │ Species12 │               0.3 │ 0.0 (0.0 - 0.0) │
    │ Species14 │               0.3 │ 0.0 (0.0 - 0.0) │
    │ Species16 │               0.3 │ 0.0 (0.0 - 0.0) │
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
    SiteB-1          │        0.0         0.0  …     0.12865         0.0
    SiteB-2          │        0.0         0.0            0.0         0.0
    SiteB-3          │        0.0         0.0       0.023333         0.0
    SiteB-4          │        0.0         0.0            0.0         0.0
    SiteB-5          │        0.0         0.0  …         0.0         0.0

Three methods will be demonstrated.

### Jaccard Similarity

### Czekanowski Index

First, let’s compose a syntopic table object from the “y” sample data
and extract the syntopic tables in matrix format.

``` julia
syn_y = VegSci.compose_syntopic_table_object("Sample", y)
syn_y_mat = VegSci.extract_syntopic_matrix(syn_y)
syn_1_mat = VegSci.extract_syntopic_matrix(syn_1)
syn_2_mat = VegSci.extract_syntopic_matrix(syn_2)
```

    1×17 Named Matrix{Float64}
    A ╲ B │   Species1    Species2  …   Species19   Species20
    ──────┼──────────────────────────────────────────────────
    Syn2  │   0.111855    0.153854  …    0.101154   0.0567083

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
     A ╲ B │   Species4    Species5  …   Species17   Species20
    ───────┼──────────────────────────────────────────────────
    Sample │   0.104638   0.0150727  …         0.0         0.0
    Syn1   │  0.0642311    0.107037      0.0981027    0.063808
    Syn2   │        0.0    0.121978  …   0.0648035   0.0567083

``` julia
VegSci.czekanowski_index(merged_syn_mats[[:"Sample"],:], merged_syn_mats[Not(:"Sample"), :])
```

    1×2 Named Matrix{Float64}
     A ╲ B │     Syn1      Syn2
    ───────┼───────────────────
    Sample │ 0.305457  0.231434

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
