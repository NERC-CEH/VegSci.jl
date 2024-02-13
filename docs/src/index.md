

``` @meta
CurrentModule = EcoVeg
```

Documentation for
[EcoVeg.jl](https://github.com/ZekeMarhshall/EcoVeg.jl).

# EcoVeg

``` julia
using Pkg; Pkg.activate("docs")
using EcoVeg
using InvertedIndices
```

      Activating project at `~/Github/EcoVeg.jl/docs/src/docs`

Tools for vegetation science.

## Background

`EcoVeg.jl` contains tools for vegetation science using the Julia
programming language (Bezanson et al. 2017).

Solves two language problem (Roesch et al. 2023)

Aims to collate functionality found in JUICE, vegan, MAVIS into a single
location with a user-friendly API and transparent methodologies. With
the aim of assisting in the creation of reproducible analysis (Sperandii
et al. 2024).

Nomenclature follows Theurillat et al. (2021).

## Installation

To install the latest stable release of `EcoVeg`:

``` julia
using Pkg
Pkg.add("EcoVeg")
```

To install the development version of `Ecoveg`:

``` julia
using Pkg
Pkg.add(url="https://github.com/ZekeMarshall/EcoVeg.jl", rev = "develop")
```

## Usage Example

To demonstrate…

First we begin with generating two example plot by species
`NamedArrays.NamedMatrix` object using the function
`EcoVeg.generate_test_array` as test data.

``` julia
x = generate_test_array(rown = 20, coln = 30, meancoloccs = 5, rowprefix = "SiteA-", colprefix = "Species")
```

    20×30 Named Matrix{Float64}
    Releve ╲ Species │   Species1    Species2  …   Species29   Species30
    ─────────────────┼──────────────────────────────────────────────────
    SiteA-1          │  0.0502419         0.0  …    0.103484         0.0
    SiteA-2          │        0.0         0.0       0.165246         0.0
    SiteA-3          │        0.0         0.0       0.369688         0.0
    SiteA-4          │        0.0         0.0            0.0         0.0
    SiteA-5          │        0.0    0.183766       0.449481         0.0
    SiteA-6          │        0.0    0.335033            0.0         0.0
    SiteA-7          │        0.0         0.0            0.0         0.0
    SiteA-8          │   0.518155         0.0            0.0         0.0
    SiteA-9          │        0.0    0.149582            0.0         0.0
    SiteA-10         │        0.0         0.0            0.0         0.0
    SiteA-11         │        0.0         0.0            0.0         0.0
    SiteA-12         │        0.0         0.0            0.0         0.0
    SiteA-13         │        0.0         0.0            0.0         0.0
    SiteA-14         │        0.0         0.0       0.304118         0.0
    SiteA-15         │   0.234075         0.0            0.0         0.0
    SiteA-16         │        1.0         0.0            0.0         0.0
    SiteA-17         │        0.0         0.0            0.0         0.0
    SiteA-18         │        0.0         0.0            0.0         0.0
    SiteA-19         │        0.0         0.0            0.0         0.0
    SiteA-20         │ 0.00645557         0.0  …   0.0266846         0.0

### Classification

Let’s artifically create some clusters for now…

``` julia
cluster1 = ["SiteA-1", "SiteA-2", "SiteA-4", "SiteA-7", "SiteA-10", "SiteA-11", "SiteA-12", "SiteA-15", "SiteA-18", "SiteA-19"]
cluster2 = ["SiteA-3", "SiteA-5", "SiteA-6", "SiteA-8", "SiteA-9", "SiteA-13", "SiteA-14", "SiteA-16", "SiteA-17", "SiteA-20"]
```

    10-element Vector{String}:
     "SiteA-3"
     "SiteA-5"
     "SiteA-6"
     "SiteA-8"
     "SiteA-9"
     "SiteA-13"
     "SiteA-14"
     "SiteA-16"
     "SiteA-17"
     "SiteA-20"

### Creation of Syntopic Tables

Once the plots have been grouped into clusters, we can proceed to
summarise their composition via the creation of `SyntopicTable` objects.

``` julia
syn_1 = EcoVeg.compose_syntopic_table_object("Syn1", x[cluster1,:])
syn_2 = EcoVeg.compose_syntopic_table_object("Syn2", x[cluster2,:])
print_summary_syntopic_table(syn_2, "normal", "cover_proportion")
```



    Community Name: Syn2
    Releves: n = 10
    Species: n = 22
    ┌───────────┬───────────────────┬─────────────────┐
    │   Species │ RelativeFrequency │       Abundance │
    ├───────────┼───────────────────┼─────────────────┤
    │ Species29 │               0.4 │ 0.3 (0.0 - 0.4) │
    │  Species1 │               0.3 │ 0.5 (0.0 - 1.0) │
    │  Species2 │               0.3 │ 0.2 (0.1 - 0.3) │
    │ Species16 │               0.2 │ 0.5 (0.1 - 1.0) │
    │ Species13 │               0.2 │ 0.4 (0.4 - 0.5) │
    │ Species27 │               0.2 │ 0.3 (0.2 - 0.5) │
    │  Species6 │               0.2 │ 0.2 (0.1 - 0.3) │
    │  Species8 │               0.2 │ 0.2 (0.0 - 0.3) │
    │ Species14 │               0.2 │ 0.2 (0.1 - 0.2) │
    │ Species23 │               0.2 │ 0.2 (0.0 - 0.4) │
    │ Species24 │               0.2 │ 0.2 (0.1 - 0.3) │
    │ Species18 │               0.2 │ 0.1 (0.0 - 0.1) │
    │ Species22 │               0.2 │ 0.1 (0.0 - 0.2) │
    │ Species21 │               0.2 │ 0.0 (0.0 - 0.1) │
    │ Species10 │               0.1 │ 0.4 (0.4 - 0.4) │
    │ Species25 │               0.1 │ 0.3 (0.3 - 0.3) │
    │  Species3 │               0.1 │ 0.2 (0.2 - 0.2) │
    │  Species7 │               0.1 │ 0.2 (0.2 - 0.2) │
    │  Species9 │               0.1 │ 0.2 (0.2 - 0.2) │
    │ Species15 │               0.1 │ 0.2 (0.2 - 0.2) │
    │ Species17 │               0.1 │ 0.1 (0.1 - 0.1) │
    │ Species12 │               0.1 │ 0.0 (0.0 - 0.0) │
    └───────────┴───────────────────┴─────────────────┘

### Identification of High-Fidelity Species

### Generation of Pseudo-Releves

### Assignment of Releves to Vegetation Classes

Let’s generate a second example matrix, consisting of sample 5 releves,
against which we want to calculate the similarity.

``` julia
y = generate_test_array(rown = 5, coln = 30, meancoloccs = 5, rowprefix = "SiteB-", colprefix = "Species")
```

    5×30 Named Matrix{Float64}
    Releve ╲ Species │    Species1     Species2  …    Species29    Species30
    ─────────────────┼──────────────────────────────────────────────────────
    SiteB-1          │         0.0          0.0  …      0.13311          0.0
    SiteB-2          │    0.093172          0.0        0.124078          0.0
    SiteB-3          │         0.0  0.000280564             0.0     0.162627
    SiteB-4          │         0.0          0.0        0.149663          0.0
    SiteB-5          │         0.0          0.0  …          0.0          0.0

Three methods will be demonstrated.

### Jaccard Similarity

### Czekanowski Index

First, let’s compose a syntopic table object from the “y” sample data
and extract the syntopic tables in matrix format.

``` julia
syn_y = EcoVeg.compose_syntopic_table_object("Sample", y)
syn_y_mat = extract_syntopic_matrix(syn_y)
syn_1_mat = extract_syntopic_matrix(syn_1)
syn_2_mat = extract_syntopic_matrix(syn_2)
```

    1×22 Named Matrix{Float64}
    A ╲ B │  Species1   Species2   Species3  …  Species25  Species27  Species29
    ──────┼────────────────────────────────────────────────────────────────────
    Syn2  │  0.518155   0.183766   0.204987  …   0.302602   0.313112   0.336903

Now we have three matrices, containg the relative frequencies of each
species present in the sample releves which constitute the
phytocoenosis’. However, each of the phytocoenosis is composed of a
different set of species, in Julia we need a helper function to merge
these matrices and ensure each matric contains each species across all
the matrices.

``` julia
merged_syn_mats = EcoVeg.merge_namedarrays([syn_y_mat, syn_1_mat, syn_2_mat])
```

    3×29 Named Matrix{Float64}
     A ╲ B │    Species1     Species2  …    Species27    Species28
    ───────┼──────────────────────────────────────────────────────
    Sample │    0.093172  0.000280564  …          0.0          0.0
    Syn1   │    0.142159          0.0         0.28964     0.326446
    Syn2   │    0.518155     0.183766  …     0.313112          0.0

``` julia
EcoVeg.czekanowski_index(merged_syn_mats[[:"Sample"],:], merged_syn_mats[Not(:"Sample"), :])
```

    1×2 Named Matrix{Float64}
     A ╲ B │     Syn1      Syn2
    ───────┼───────────────────
    Sample │ 0.371491  0.373745

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
