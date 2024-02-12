

``` @meta
CurrentModule = EcoVeg
```

Documentation for
[EcoVeg.jl](https://github.com/ZekeMarhshall/EcoVeg.jl).

# EcoVeg

Tools for vegetation science.

## Background

`EcoVeg.jl` contains tools for vegetation science using the Julia
programming language \[@bezanson2017\].

Solves two language problem \[@roesch2023\]

Aims to collate functionality found in JUICE, vegan, MAVIS into a single
location with a user-friendly API and transparent methodologies. With
the aim of assisting in the creation of reproducible analysis
\[@sperandii2024\].

Nomenclature follows @theurillat2021.

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
x = generate_test_array(rown = 20, coln = 30, meancoloccs = 10, rowprefix = "SiteA-", colprefix = "Species")
```

    20×30 Named Matrix{Float64}
    Releve ╲ Species │   Species1    Species2  …   Species29   Species30
    ─────────────────┼──────────────────────────────────────────────────
    SiteA-1          │        0.0         0.0  …         0.0         0.0
    SiteA-2          │        0.0         0.0            0.0   0.0742907
    SiteA-3          │  0.0704807   0.0369201            0.0         0.0
    SiteA-4          │        0.0         0.0            0.0         0.0
    SiteA-5          │        0.0         0.0       0.065016    0.114198
    SiteA-6          │        0.0         0.0     0.00555444         0.0
    SiteA-7          │        0.0   0.0550128            0.0   0.0563497
    SiteA-8          │        0.0         0.0            0.0   0.0475957
    SiteA-9          │        0.0   0.0139922            0.0         0.0
    SiteA-10         │   0.048195    0.058435      0.0296611         0.0
    SiteA-11         │        0.0         0.0            0.0         0.0
    SiteA-12         │        0.0         0.0            0.0   0.0508891
    SiteA-13         │        0.0   0.0798633            0.0    0.131321
    SiteA-14         │        0.0         0.0            0.0    0.183577
    SiteA-15         │        0.0         0.0       0.140305   0.0119757
    SiteA-16         │  0.0752609    0.067526      0.0843943         0.0
    SiteA-17         │  0.0864513         0.0            0.0         0.0
    SiteA-18         │ 0.00308486         0.0            0.0    0.147017
    SiteA-19         │        0.0    0.133836            0.0         0.0
    SiteA-20         │        0.0         0.0  …         0.0         0.0

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
print_summary_syntopic_table(syn_2)
```



    Community Name: Syn2
    Releves: n = 10
    Species: n = 30
    ╭───────────┬───────────────────┬─────────────────╮
    │   Species │ RelativeFrequency │       Abundance │
    ├───────────┼───────────────────┼─────────────────┤
    │  Species8 │          0.116874 │ 0.1 (0.1 - 0.3) │
    │ Species15 │         0.0623607 │ 0.1 (0.0 - 0.2) │
    │ Species19 │         0.0584902 │ 0.1 (0.1 - 0.2) │
    │ Species11 │         0.0523788 │ 0.1 (0.1 - 0.2) │
    │ Species28 │         0.0513153 │ 0.2 (0.0 - 0.2) │
    │ Species26 │         0.0478617 │ 0.1 (0.1 - 0.2) │
    │ Species30 │         0.0476692 │ 0.1 (0.0 - 0.2) │
    │ Species14 │         0.0445658 │ 0.1 (0.0 - 0.2) │
    │ Species18 │         0.0403507 │ 0.1 (0.1 - 0.3) │
    │ Species21 │         0.0364615 │ 0.1 (0.1 - 0.2) │
    │ Species12 │         0.0343615 │ 0.1 (0.0 - 0.2) │
    │ Species17 │         0.0325703 │ 0.1 (0.0 - 0.1) │
    │ Species16 │         0.0315202 │ 0.1 (0.0 - 0.1) │
    │ Species23 │         0.0300857 │ 0.1 (0.0 - 0.1) │
    │  Species4 │         0.0298907 │ 0.1 (0.1 - 0.2) │
    │ Species13 │         0.0284227 │ 0.1 (0.0 - 0.1) │
    │  Species5 │         0.0275004 │ 0.1 (0.0 - 0.1) │
    │ Species22 │         0.0267309 │ 0.1 (0.0 - 0.1) │
    │  Species1 │         0.0232193 │ 0.1 (0.1 - 0.1) │
    │ Species27 │         0.0223382 │ 0.1 (0.1 - 0.1) │
    │  Species2 │         0.0198302 │ 0.1 (0.0 - 0.1) │
    │  Species7 │         0.0198162 │ 0.1 (0.1 - 0.1) │
    │  Species6 │          0.017789 │ 0.1 (0.1 - 0.1) │
    │ Species20 │         0.0175654 │ 0.1 (0.1 - 0.1) │
    │ Species25 │         0.0160109 │ 0.1 (0.0 - 0.1) │
    │  Species9 │         0.0156659 │ 0.1 (0.0 - 0.1) │
    │ Species29 │         0.0154965 │ 0.1 (0.0 - 0.1) │
    │ Species10 │         0.0115231 │ 0.1 (0.0 - 0.1) │
    │  Species3 │         0.0113954 │ 0.1 (0.0 - 0.1) │
    │ Species24 │        0.00993922 │ 0.0 (0.0 - 0.1) │
    ╰───────────┴───────────────────┴─────────────────╯

### Identification of High-Fidelity Species

### Generation of Pseudo-Releves

### Assignment of Releves to Vegetation Classes

Let’s generate a second example matrix, consisting of sample 5 releves,
against which we want to calculate the similarity.

``` julia
y = generate_test_array(rown = 5, coln = 30, meancoloccs = 10, rowprefix = "SiteB-", colprefix = "Species")
```

    5×30 Named Matrix{Float64}
    Releve ╲ Species │    Species1     Species2  …    Species29    Species30
    ─────────────────┼──────────────────────────────────────────────────────
    SiteB-1          │         0.0     0.267636  …          0.0          0.0
    SiteB-2          │    0.149368    0.0515994             0.0     0.135862
    SiteB-3          │         0.0    0.0404588             0.0     0.108958
    SiteB-4          │         0.0          0.0        0.012553    0.0916869
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

    1×30 Named Matrix{Float64}
    A ╲ B │   Species8   Species15  …    Species3   Species24
    ──────┼──────────────────────────────────────────────────
    Syn2  │   0.116874   0.0623607  …   0.0113954  0.00993922

Now we have three matrices, containg the relative frequencies of each
species present in the sample releves which constitute the
phytocoenosis’. However, each of the phytocoenosis is composed of a
different set of species, in Julia we need a helper function to ensure
all columns are present in each of the matrices before joining.

``` julia
# align_array_columns(syn_y_mat, syn_2_mat)
# align_array_columns(syn_y_mat, syn_1_mat)
```

### Multivariate Analysis

### Ecological Trajectory Analysis

## External Resources

## Implemented Methodologies

## Contribute

## Acknowledgements

## References
