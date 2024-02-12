

\[![version](https://img.shields.io/badge/version-0.1-blue)\]
[![CI](https://github.com/JuliaCI/PkgTemplates.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ZekeMarshall/EcoVeg.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Codecov](https://codecov.io/gh/ZekeMarshall/EcoVeg.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ZekeMarshall/EcoVeg.jl)
[![Aqua
QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/ZekeMarshall/EcoVeg.jl)
[![License](https://img.shields.io/github/license/ZekeMarshall/EcoVeg.jl)](LICENSE)

# EcoVeg

Tools for vegetation science.

## Background

in vegetation science using the Julia programming language (Bezanson et
al. 2017)

Solves two language problem (Roesch et al. 2023)

Aims to collate functionality found in JUICE, vegan, MAVIS into a single
location with a user-friendly API and transparent methodologies.

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
x = generate_test_array(rown = 20, coln = 30, zerop = 0.6, rowprefix = "SiteA-", colprefix = "Species")
```

### Classification

Let’s artifically create some clusters for now…

``` julia
cluster1 = ["SiteA-1", "SiteA-2", "SiteA-4", "SiteA-7", "SiteA-10", "SiteA-11", "SiteA-12", "SiteA-15", "SiteA-18", "SiteA-19"]
cluster2 = ["SiteA-3", "SiteA-5", "SiteA-6", "SiteA-8", "SiteA-9", "SiteA-13", "SiteA-14", "SiteA-16", "SiteA-17", "SiteA-20"]
```

### Creation of Syntopic Tables

Once the plots have been grouped into clusters, we can proceed to
summarise their composition via the creation of `SyntopicTable` objects.

``` julia
syn_1 = EcoVeg.compose_syntopic_table_object("Syn1", x[cluster1,:])
syn_2 = EcoVeg.compose_syntopic_table_object("Syn2", x[cluster2,:])
print_summary_syntopic_table(syn_2)
```

### Identification of High-Fidelity Species

### Generation of Pseudo-Releves

### Assignment of Releves to Vegetation Classes

Let’s generate a second example matrix, consisting of sample 5 releves,
against which we want to calculate the similarity.

``` julia
y = generate_test_array(rown = 5, coln = 30, zerop = 0.6, rowprefix = "SiteB-", colprefix = "Species")
```

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

Theurillat, Jean-Paul, Wolfgang Willner, Federico Fernández-González,
Helga Bültmann, Andraž Čarni, Daniela Gigante, Ladislav Mucina, and
Heinrich Weber. 2021. “International Code of Phytosociological
Nomenclature. 4th Edition.” *Applied Vegetation Science* 24 (1): e12491.
<https://doi.org/10.1111/avsc.12491>.
