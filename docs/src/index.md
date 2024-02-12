

``` @meta
CurrentModule = EcoVeg
```

Documentation for
[EcoVeg.jl](https://github.com/ZekeMarhshall/EcoVeg.jl).

# EcoVeg

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
x = generate_test_array(rown = 20, coln = 30, meancoloccs = 10, rowprefix = "SiteA-", colprefix = "Species")
```

    20×30 Named Matrix{Float64}
    Releve ╲ Species │    Species1     Species2  …    Species29    Species30
    ─────────────────┼──────────────────────────────────────────────────────
    SiteA-1          │         0.0          0.0  …          0.0          0.0
    SiteA-2          │         0.0          0.0       0.0782331          0.0
    SiteA-3          │         0.0    0.0575326        0.155053     0.152938
    SiteA-4          │         0.0     0.144142        0.128875          0.0
    SiteA-5          │         0.0          0.0             0.0          0.0
    SiteA-6          │    0.085319          0.0             0.0          0.0
    SiteA-7          │         0.0          0.0             0.0          0.0
    SiteA-8          │   0.0449193          0.0             0.0          0.0
    SiteA-9          │         0.0          0.0        0.216531          0.0
    SiteA-10         │         0.0          0.0             0.0      0.18558
    SiteA-11         │   0.0771162          0.0        0.242951          0.0
    SiteA-12         │         0.0          0.0             0.0          0.0
    SiteA-13         │         0.0          0.0             0.0     0.329161
    SiteA-14         │         0.0     0.177045             0.0      0.14028
    SiteA-15         │    0.120292          0.0             0.0      0.28692
    SiteA-16         │   0.0241883     0.111589       0.0757173          0.0
    SiteA-17         │    0.138463          0.0             0.0     0.155279
    SiteA-18         │         0.0          0.0             0.0     0.204505
    SiteA-19         │    0.114521    0.0774826             0.0          0.0
    SiteA-20         │    0.108287          0.0  …    0.0175558    0.0594922

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
    │ Species14 │         0.0974639 │ 0.1 (0.1 - 0.3) │
    │ Species30 │          0.083715 │ 0.2 (0.1 - 0.3) │
    │ Species20 │         0.0689173 │ 0.1 (0.0 - 0.2) │
    │  Species8 │         0.0600403 │ 0.1 (0.0 - 0.3) │
    │ Species19 │         0.0497797 │ 0.1 (0.1 - 0.2) │
    │ Species27 │         0.0493997 │ 0.2 (0.1 - 0.2) │
    │ Species29 │         0.0464857 │ 0.1 (0.0 - 0.2) │
    │ Species23 │         0.0424856 │ 0.1 (0.1 - 0.2) │
    │  Species1 │         0.0401177 │ 0.1 (0.0 - 0.1) │
    │ Species16 │         0.0376498 │ 0.1 (0.0 - 0.2) │
    │ Species22 │         0.0363398 │ 0.1 (0.0 - 0.2) │
    │ Species12 │         0.0353443 │ 0.2 (0.2 - 0.2) │
    │  Species2 │         0.0346166 │ 0.1 (0.1 - 0.2) │
    │ Species28 │         0.0339354 │ 0.1 (0.0 - 0.2) │
    │ Species17 │         0.0302904 │ 0.1 (0.0 - 0.1) │
    │  Species4 │          0.027939 │ 0.1 (0.0 - 0.2) │
    │ Species21 │         0.0260201 │ 0.1 (0.1 - 0.1) │
    │ Species13 │         0.0255388 │ 0.0 (0.0 - 0.1) │
    │  Species3 │         0.0246532 │ 0.0 (0.0 - 0.1) │
    │ Species15 │         0.0245426 │ 0.1 (0.0 - 0.1) │
    │ Species26 │         0.0197612 │ 0.0 (0.0 - 0.1) │
    │ Species11 │         0.0182556 │ 0.1 (0.1 - 0.1) │
    │ Species24 │         0.0169718 │ 0.1 (0.1 - 0.1) │
    │ Species10 │         0.0155926 │ 0.1 (0.0 - 0.1) │
    │  Species7 │         0.0151053 │ 0.1 (0.1 - 0.1) │
    │  Species5 │         0.0123686 │ 0.1 (0.0 - 0.1) │
    │ Species18 │         0.0112344 │ 0.1 (0.1 - 0.1) │
    │  Species9 │        0.00908157 │ 0.1 (0.1 - 0.1) │
    │  Species6 │        0.00406966 │ 0.0 (0.0 - 0.0) │
    │ Species25 │        0.00228414 │ 0.0 (0.0 - 0.0) │
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
    Releve ╲ Species │   Species1    Species2  …   Species29   Species30
    ─────────────────┼──────────────────────────────────────────────────
    SiteB-1          │        0.0         0.0  …         0.0         0.0
    SiteB-2          │        0.0   0.0307212       0.232576         0.0
    SiteB-3          │        0.0         0.0       0.046245   0.0969491
    SiteB-4          │   0.081061         0.0      0.0764102         0.0
    SiteB-5          │   0.101557         0.0  …         0.0         0.0

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
    A ╲ B │  Species14   Species30  …    Species6   Species25
    ──────┼──────────────────────────────────────────────────
    Syn2  │  0.0974639    0.083715  …  0.00406966  0.00228414

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

Bezanson, Jeff, Alan Edelman, Stefan Karpinski, and Viral B Shah. 2017.
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
