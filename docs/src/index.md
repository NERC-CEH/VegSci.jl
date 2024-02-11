

``` @meta
CurrentModule = EcoVeg
```

<!-- ![](assets/wide_logo.png) -->

Documentation for
[EcoVeg.jl](https://github.com/ZekeMarhshall/EcoVeg.jl).

# EcoVeg

      Activating new project at `~/Github/EcoVeg.jl/docs/src/docs`

`EcoVeg.jl` is a package containing functions for the analysis of
vegetation plot sample data….

in vegetation science using the Julia programming language
\[@bezanson2017\]

Solves two language problem \[@roesch2023\]

Aims to collate functionality found in JUICE, vegan, MAVIS into a single
location with a user-friendly API and transparent methodologies.

## Example Pluto Notebook

### Local Tour

To run the tour locally, just clone this repo and start `Pluto.jl` as
follows:

``` julia
] add Pluto
using Pluto
Pluto.run()
```

All notebooks are contained in `docs/pluto`.

## Background

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
x = generate_test_array(rown = 20, coln = 30, zerop = 0.6, rowprefix = "SiteA-", colprefix = "Species")
```

    20×30 Named Matrix{Float64}
    Releve ╲ Species │    Species1     Species2  …    Species29    Species30
    ─────────────────┼──────────────────────────────────────────────────────
    SiteA-1          │    0.574341     0.835083  …          0.0     0.647713
    SiteA-2          │    0.308076     0.952746        0.491991          0.0
    SiteA-3          │   0.0170512          0.0             0.0     0.995746
    SiteA-4          │         0.0     0.182911         0.98594     0.902697
    SiteA-5          │  0.00449277     0.978113        0.950479      0.19102
    SiteA-6          │         0.0          0.0        0.103873     0.780704
    SiteA-7          │         0.0     0.807792        0.245446          0.0
    SiteA-8          │   0.0865799     0.785077             0.0          0.0
    SiteA-9          │    0.592954     0.443407             0.0          0.0
    SiteA-10         │    0.940073     0.766005       0.0214637     0.677974
    SiteA-11         │    0.625696          0.0        0.376025     0.655731
    SiteA-12         │    0.279377     0.369117        0.768802     0.851647
    SiteA-13         │    0.352045     0.218777             0.0          0.0
    SiteA-14         │         0.0     0.497137        0.476686     0.594838
    SiteA-15         │         0.0          0.0             0.0     0.114054
    SiteA-16         │    0.774778     0.538604        0.910406          0.0
    SiteA-17         │         0.0     0.811507             0.0     0.701075
    SiteA-18         │    0.639591          0.0             0.0          0.0
    SiteA-19         │     0.55452     0.860287             0.0          0.0
    SiteA-20         │     0.14533     0.595273  …     0.523479     0.907187

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
    │ Species19 │          0.523543 │ 0.8 (0.2 - 0.9) │
    │  Species2 │          0.486789 │ 0.6 (0.2 - 1.0) │
    │  Species4 │          0.462055 │ 0.9 (0.4 - 1.0) │
    │ Species15 │           0.44981 │ 0.7 (0.0 - 1.0) │
    │ Species27 │          0.434385 │ 0.6 (0.2 - 1.0) │
    │ Species30 │          0.417057 │ 0.7 (0.2 - 1.0) │
    │  Species6 │          0.405933 │ 0.5 (0.0 - 0.8) │
    │ Species24 │          0.392898 │ 0.5 (0.1 - 0.9) │
    │ Species10 │          0.381608 │ 0.4 (0.1 - 0.9) │
    │ Species21 │          0.378836 │ 0.7 (0.2 - 1.0) │
    │ Species14 │          0.369817 │ 0.5 (0.1 - 1.0) │
    │ Species28 │          0.353245 │ 0.5 (0.2 - 1.0) │
    │ Species22 │          0.336727 │ 0.5 (0.2 - 1.0) │
    │ Species20 │          0.322975 │ 0.7 (0.2 - 0.9) │
    │  Species5 │          0.314488 │ 0.7 (0.2 - 1.0) │
    │ Species18 │          0.307546 │ 0.5 (0.0 - 0.9) │
    │  Species8 │          0.303783 │ 0.4 (0.0 - 0.8) │
    │ Species29 │          0.296492 │ 0.5 (0.1 - 1.0) │
    │  Species3 │          0.293004 │ 0.8 (0.5 - 0.8) │
    │ Species25 │          0.284239 │ 0.6 (0.1 - 1.0) │
    │ Species26 │          0.281381 │ 0.5 (0.1 - 0.7) │
    │  Species7 │           0.26627 │ 0.3 (0.0 - 0.8) │
    │ Species16 │          0.261026 │ 0.5 (0.1 - 0.8) │
    │ Species12 │           0.25653 │ 0.4 (0.1 - 1.0) │
    │ Species23 │          0.224693 │ 0.7 (0.1 - 0.8) │
    │ Species17 │          0.221477 │ 0.4 (0.1 - 0.9) │
    │  Species1 │          0.197323 │ 0.1 (0.0 - 0.8) │
    │  Species9 │          0.178729 │ 0.2 (0.1 - 0.8) │
    │ Species11 │          0.136716 │ 0.7 (0.4 - 1.0) │
    │ Species13 │          0.126797 │ 0.2 (0.0 - 0.7) │
    ╰───────────┴───────────────────┴─────────────────╯

### Identification of High-Fidelity Species

### Generation of Pseudo-Releves

### Assignment of Releves to Vegetation Classes

Let’s generate a second example matrix, consisting of sample 5 releves,
against which we want to calculate the similarity.

``` julia
y = generate_test_array(rown = 5, coln = 30, zerop = 0.6, rowprefix = "SiteB-", colprefix = "Species")
```

    5×30 Named Matrix{Float64}
    Releve ╲ Species │  Species1   Species2  …  Species29  Species30
    ─────────────────┼──────────────────────────────────────────────
    SiteB-1          │  0.270629   0.888484  …  0.0570284  0.0940625
    SiteB-2          │       0.0   0.127423      0.303841        0.0
    SiteB-3          │       0.0   0.944684     0.0229587   0.756141
    SiteB-4          │   0.32658        0.0      0.439616   0.486081
    SiteB-5          │  0.978606        0.0  …   0.399057   0.139963

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
    A ╲ B │ Species19   Species2   Species4  …   Species9  Species11  Species13
    ──────┼────────────────────────────────────────────────────────────────────
    Syn2  │  0.523543   0.486789   0.462055  …   0.178729   0.136716   0.126797

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
