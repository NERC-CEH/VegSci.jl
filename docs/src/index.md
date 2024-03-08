

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
    SiteA-1          │         0.0          0.0  …     0.157852   0.00555798
    SiteA-2          │    0.130109          0.0       0.0747752          0.0
    SiteA-3          │    0.153157     0.177941             0.0    0.0164033
    SiteA-4          │         0.0     0.108419        0.102789    0.0942697
    SiteA-5          │   0.0104223   0.00789441        0.152983     0.139453
    SiteA-6          │         0.0    0.0747113       0.0036627     0.154264
    SiteA-7          │         0.0     0.104619        0.121775    0.0622911
    SiteA-8          │         0.0          0.0        0.100688     0.054472
    SiteA-9          │         0.0          0.0             0.0          0.0
    SiteA-10         │    0.180499          0.0             0.0          0.0
    SiteA-11         │    0.164635          0.0        0.161734     0.187418
    ⋮                            ⋮            ⋮  ⋱            ⋮            ⋮
    SiteA-20         │   0.0763144          0.0        0.105427          0.0
    SiteA-21         │    0.015283    0.0269463       0.0264219          0.0
    SiteA-22         │         0.0          0.0             0.0     0.124375
    SiteA-23         │    0.127401          0.0        0.177078          0.0
    SiteA-24         │         0.0      0.16006        0.192652     0.112483
    SiteA-25         │   0.0192523     0.114235       0.0634028     0.119082
    SiteA-26         │   0.0798989          0.0       0.0466872          0.0
    SiteA-27         │   0.0228354          0.0       0.0488663     0.103179
    SiteA-28         │         0.0    0.0207325             0.0          0.0
    SiteA-29         │         0.0          0.0             0.0          0.0
    SiteA-30         │         0.0          0.0  …          0.0          0.0

### Classification

Let’s identify some clusters, storing the cluster-releve memberships as
a dictionary.

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
      "1" => ["SiteA-4", "SiteA-5", "SiteA-6", "SiteA-11", "SiteA-12", "SiteA-14", …
      "2" => ["SiteA-2", "SiteA-3", "SiteA-7", "SiteA-15", "SiteA-21", "SiteA-25"]
      "3" => ["SiteA-1", "SiteA-8", "SiteA-9", "SiteA-10", "SiteA-13", "SiteA-17", …

### Creation of Syntopic Tables

Once the plots have been grouped into clusters, we can proceed to
summarise their composition via the creation of `SyntopicTable` objects.

``` julia
syn_1 = VegSci.compose_syntopic_table_object("Syn1", x[getindex(clusters, "1"),:])
syn_2 = VegSci.compose_syntopic_table_object("Syn2", x[getindex(clusters, "2"),:])
VegSci.print_summary_syntopic_table(syn_2, "normal", "cover_proportion")
```



    Community Name: Syn2
    Releves: n = 6
    Species: n = 19
    ┌───────────┬───────────────────┬─────────────────┐
    │   Species │ RelativeFrequency │       Abundance │
    ├───────────┼───────────────────┼─────────────────┤
    │  Species3 │               1.0 │ 0.1 (0.0 - 0.2) │
    │ Species16 │               1.0 │ 0.1 (0.1 - 0.1) │
    │  Species1 │               0.8 │ 0.1 (0.0 - 0.2) │
    │  Species2 │               0.8 │ 0.1 (0.0 - 0.2) │
    │  Species7 │               0.8 │ 0.1 (0.0 - 0.1) │
    │ Species12 │               0.8 │ 0.1 (0.0 - 0.1) │
    │  Species4 │               0.7 │ 0.1 (0.0 - 0.1) │
    │  Species5 │               0.7 │ 0.1 (0.0 - 0.1) │
    │ Species10 │               0.7 │ 0.1 (0.1 - 0.1) │
    │ Species15 │               0.7 │ 0.1 (0.0 - 0.1) │
    │ Species18 │               0.7 │ 0.1 (0.1 - 0.1) │
    │ Species19 │               0.7 │ 0.1 (0.0 - 0.1) │
    │ Species11 │               0.5 │ 0.1 (0.0 - 0.1) │
    │ Species17 │               0.5 │ 0.1 (0.0 - 0.1) │
    │ Species20 │               0.5 │ 0.1 (0.0 - 0.1) │
    │ Species14 │               0.3 │ 0.0 (0.0 - 0.0) │
    │  Species6 │               0.2 │ 0.1 (0.1 - 0.1) │
    │  Species8 │               0.2 │ 0.0 (0.0 - 0.0) │
    │  Species9 │               0.2 │ 0.0 (0.0 - 0.0) │
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
    SiteB-1          │        0.0         0.0  …         0.0         0.0
    SiteB-2          │        0.0   0.0641312      0.0955262         0.0
    SiteB-3          │        0.0         0.0       0.045238         0.0
    SiteB-4          │        0.0         0.0            0.0         0.0
    SiteB-5          │        0.0    0.245125  …    0.181263         0.0

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

    1×19 Named Matrix{Float64}
    A ╲ B │   Species1    Species2  …   Species19   Species20
    ──────┼──────────────────────────────────────────────────
    Syn2  │   0.130109    0.114235  …    0.069089   0.0622911

Now we have three matrices, containg the relative frequencies of each
species present in the sample releves which constitute each syntaxon.
However, each of the syntaxa are composed of a different set of species,
in Julia we need a helper function to merge these matrices and ensure
each matrix contains each species across all the matrices. This function
is broadly equivalent to the R function `base::merge`.

``` julia
merged_syn_mats = VegSci.merge_namedarrays([syn_y_mat, syn_1_mat, syn_2_mat])
```

    3×26 Named Matrix{Float64}
     A ╲ B │   Species2    Species3  …   Species19   Species20
    ───────┼──────────────────────────────────────────────────
    Sample │   0.154628    0.130323  …         0.0         0.0
    Syn1   │  0.0860061   0.0715421       0.104108    0.134974
    Syn2   │   0.114235    0.119714  …    0.069089   0.0622911

``` julia
VegSci.steinhaus_coefficient(merged_syn_mats[[:"Sample"],:], merged_syn_mats[Not(:"Sample"), :])
```

    1×2 Named Matrix{Float64}
     A ╲ B │     Syn1      Syn2
    ───────┼───────────────────
    Sample │ 0.302282  0.268195

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
