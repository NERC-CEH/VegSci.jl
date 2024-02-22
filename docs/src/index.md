

``` @meta
CurrentModule = VegSci
```

Documentation for
[VegSci.jl](https://github.com/ZekeMarhshall/VegSci.jl).

# VegSci

Tools for vegetation science.

## Background

`VegSci.jl` contains tools for vegetation science using the Julia
programming language (Bezanson et al. 2017).

Solves two language problem (Roesch et al. 2023)

Aims to collate functionality found in JUICE, vegan, MAVIS into a single
location with a user-friendly API and transparent methodologies. With
the aim of assisting in the creation of reproducible analysis (Sperandii
et al. 2024).

Nomenclature follows Theurillat et al. (2021).

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
    SiteA-1          │   0.126303    0.187072  …    0.210793         0.0
    SiteA-2          │        0.0    0.144847            0.0   0.0936073
    SiteA-3          │   0.125143   0.0774974            0.0   0.0216059
    SiteA-4          │        0.0         0.0            0.0   0.0820958
    SiteA-5          │    0.11744         0.0            0.0         0.0
    SiteA-6          │   0.143738         0.0            0.0   0.0366179
    SiteA-7          │   0.200434         0.0       0.109825         0.0
    SiteA-8          │   0.039486   0.0575995            0.0   0.0494997
    SiteA-9          │        0.0         0.0        0.16484   0.0646168
    SiteA-10         │   0.044275   0.0280082            0.0         0.0
    SiteA-11         │   0.112592    0.144315            0.0         0.0
    ⋮                           ⋮           ⋮  ⋱           ⋮           ⋮
    SiteA-20         │   0.257369    0.155278            0.0         0.0
    SiteA-21         │   0.187503   0.0431709            0.0   0.0543969
    SiteA-22         │        0.0         0.0       0.183636         0.0
    SiteA-23         │        0.0    0.139459       0.094891         0.0
    SiteA-24         │  0.0755919         0.0       0.042628    0.143788
    SiteA-25         │        0.0    0.147187            0.0         0.0
    SiteA-26         │  0.0264285         0.0       0.152865         0.0
    SiteA-27         │        0.0    0.246575       0.209535    0.188035
    SiteA-28         │        0.0         0.0            0.0         0.0
    SiteA-29         │   0.117252         0.0      0.0344899         0.0
    SiteA-30         │        0.0   0.0408142  …         0.0   0.0131516

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
      2 => [1, 2, 5, 7, 17, 18, 19, 21, 23, 26, 27, 29, 30]
      3 => [3, 10, 11, 13, 20, 25]
      1 => [4, 6, 8, 9, 12, 14, 15, 16, 22, 24, 28]

### Creation of Syntopic Tables

Once the plots have been grouped into clusters, we can proceed to
summarise their composition via the creation of `SyntopicTable` objects.

``` julia
syn_1 = VegSci.compose_syntopic_table_object("Syn1", x[getindex(memberships, 1),:])
syn_2 = VegSci.compose_syntopic_table_object("Syn2", x[getindex(memberships, 2),:])
VegSci.print_summary_syntopic_table(syn_2, "normal", "cover_proportion")
```



    Community Name: Syn2
    Releves: n = 13
    Species: n = 20
    ┌───────────┬───────────────────┬─────────────────┐
    │   Species │ RelativeFrequency │       Abundance │
    ├───────────┼───────────────────┼─────────────────┤
    │ Species16 │               0.8 │ 0.1 (0.0 - 0.2) │
    │ Species15 │               0.7 │ 0.1 (0.1 - 0.3) │
    │  Species1 │               0.6 │ 0.1 (0.0 - 0.3) │
    │  Species2 │               0.6 │ 0.1 (0.0 - 0.2) │
    │  Species4 │               0.6 │ 0.1 (0.0 - 0.2) │
    │ Species12 │               0.6 │ 0.1 (0.0 - 0.2) │
    │ Species14 │               0.6 │ 0.1 (0.0 - 0.1) │
    │ Species19 │               0.6 │ 0.1 (0.0 - 0.2) │
    │  Species3 │               0.5 │ 0.2 (0.1 - 0.2) │
    │  Species5 │               0.5 │ 0.1 (0.0 - 0.2) │
    │  Species8 │               0.5 │ 0.0 (0.0 - 0.1) │
    │ Species18 │               0.5 │ 0.0 (0.0 - 0.1) │
    │ Species10 │               0.4 │ 0.1 (0.0 - 0.2) │
    │ Species17 │               0.4 │ 0.1 (0.0 - 0.1) │
    │ Species20 │               0.4 │ 0.1 (0.0 - 0.2) │
    │  Species6 │               0.3 │ 0.1 (0.0 - 0.1) │
    │  Species7 │               0.2 │ 0.1 (0.1 - 0.2) │
    │ Species11 │               0.2 │ 0.1 (0.1 - 0.1) │
    │ Species13 │               0.2 │ 0.1 (0.0 - 0.2) │
    │  Species9 │               0.2 │ 0.0 (0.0 - 0.1) │
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
    SiteB-1          │        0.0         0.0  …    0.457253         0.0
    SiteB-2          │   0.393232         0.0            0.0         0.0
    SiteB-3          │   0.575143         0.0       0.055049         0.0
    SiteB-4          │        0.0         0.0            0.0         0.0
    SiteB-5          │        0.0         0.0  …         0.0   0.0581744

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

    1×20 Named Matrix{Float64}
    A ╲ B │  Species1   Species2   Species3  …  Species18  Species19  Species20
    ──────┼────────────────────────────────────────────────────────────────────
    Syn2  │  0.121872   0.135754   0.174128  …  0.0348915   0.131345  0.0619423

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
     A ╲ B │  Species1   Species3   Species5  …  Species14  Species15  Species19
    ───────┼────────────────────────────────────────────────────────────────────
    Sample │  0.484188   0.185498   0.242788  …        0.0        0.0        0.0
    Syn1   │ 0.0554367  0.0892419   0.094885     0.0875246  0.0337278   0.107604
    Syn2   │  0.121872   0.174128   0.132866  …  0.0861042   0.142692   0.131345

``` julia
VegSci.czekanowski_index(merged_syn_mats[[:"Sample"],:], merged_syn_mats[Not(:"Sample"), :])
```

    1×2 Named Matrix{Float64}
     A ╲ B │     Syn1      Syn2
    ───────┼───────────────────
    Sample │ 0.393623  0.467868

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
