

``` @meta
CurrentModule = VegSci
```

Documentation for
[VegSci.jl](https://github.com/ZekeMarhshall/VegSci.jl).

# VegSci

``` julia
using VegSci
using InvertedIndices
using Clustering
```

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
x = generate_test_array(rown = 30, coln = 20, meancoloccs = 10, rowprefix = "SiteA-", colprefix = "Species")
```

    30×20 Named Matrix{Float64}
    Releve ╲ Species │    Species1     Species2  …    Species19    Species20
    ─────────────────┼──────────────────────────────────────────────────────
    SiteA-1          │         0.0     0.178572  …          0.0          0.0
    SiteA-2          │         0.0          0.0       0.0656326     0.150752
    SiteA-3          │    0.127538          0.0             0.0          0.0
    SiteA-4          │   0.0865866    0.0525902       0.0265362          0.0
    SiteA-5          │   0.0863782     0.153892      0.00376953          0.0
    SiteA-6          │    0.122489    0.0869498        0.127231          0.0
    SiteA-7          │         0.0     0.228873             0.0    0.0726601
    SiteA-8          │   0.0261627  0.000347997       0.0771068    0.0120358
    SiteA-9          │         0.0          0.0       0.0205265    0.0539659
    SiteA-10         │    0.106896    0.0260744             0.0          0.0
    SiteA-11         │         0.0          0.0       0.0389704          0.0
    ⋮                            ⋮            ⋮  ⋱            ⋮            ⋮
    SiteA-20         │   0.0889668     0.160342             0.0     0.147702
    SiteA-21         │         0.0     0.219887        0.108587    0.0439591
    SiteA-22         │    0.282867          0.0             0.0    0.0424041
    SiteA-23         │   0.0277427          0.0       0.0763637    0.0939302
    SiteA-24         │         0.0     0.177612             0.0          0.0
    SiteA-25         │         0.0          0.0             0.0          0.0
    SiteA-26         │         0.0     0.167282        0.162137          0.0
    SiteA-27         │    0.235794          0.0       0.0409572          0.0
    SiteA-28         │         0.0     0.111215       0.0519763     0.131982
    SiteA-29         │   0.0479258          0.0             0.0          0.0
    SiteA-30         │   0.0444138    0.0955813  …     0.115388    0.0969775

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
      2 => [5, 16, 24, 25]
      3 => [2, 6, 7, 9, 14, 17, 18, 20, 21, 26, 28, 30]
      1 => [1, 3, 4, 8, 10, 11, 12, 13, 15, 19, 22, 23, 27, 29]

### Creation of Syntopic Tables

Once the plots have been grouped into clusters, we can proceed to
summarise their composition via the creation of `SyntopicTable` objects.

``` julia
syn_1 = VegSci.compose_syntopic_table_object("Syn1", x[getindex(memberships, 1),:])
syn_2 = VegSci.compose_syntopic_table_object("Syn2", x[getindex(memberships, 2),:])
print_summary_syntopic_table(syn_2, "normal", "cover_proportion")
```



    Community Name: Syn2
    Releves: n = 4
    Species: n = 18
    ┌───────────┬───────────────────┬─────────────────┐
    │   Species │ RelativeFrequency │       Abundance │
    ├───────────┼───────────────────┼─────────────────┤
    │ Species11 │               0.8 │ 0.2 (0.2 - 0.3) │
    │  Species4 │               0.8 │ 0.1 (0.0 - 0.2) │
    │ Species16 │               0.8 │ 0.1 (0.1 - 0.3) │
    │  Species2 │               0.5 │ 0.2 (0.2 - 0.2) │
    │  Species3 │               0.5 │ 0.2 (0.2 - 0.2) │
    │  Species8 │               0.5 │ 0.2 (0.1 - 0.2) │
    │  Species9 │               0.5 │ 0.2 (0.1 - 0.3) │
    │ Species14 │               0.5 │ 0.1 (0.1 - 0.1) │
    │ Species18 │               0.5 │ 0.1 (0.0 - 0.1) │
    │ Species17 │               0.5 │ 0.0 (0.0 - 0.0) │
    │ Species10 │               0.2 │ 0.2 (0.2 - 0.2) │
    │ Species12 │               0.2 │ 0.2 (0.2 - 0.2) │
    │ Species20 │               0.2 │ 0.2 (0.2 - 0.2) │
    │  Species1 │               0.2 │ 0.1 (0.1 - 0.1) │
    │ Species15 │               0.2 │ 0.1 (0.1 - 0.1) │
    │  Species6 │               0.2 │ 0.0 (0.0 - 0.0) │
    │ Species13 │               0.2 │ 0.0 (0.0 - 0.0) │
    │ Species19 │               0.2 │ 0.0 (0.0 - 0.0) │
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
    Releve ╲ Species │   Species1    Species2  …   Species29   Species30
    ─────────────────┼──────────────────────────────────────────────────
    SiteB-1          │        0.0         0.0  …         0.0         0.0
    SiteB-2          │        0.0         0.0            0.0         1.0
    SiteB-3          │        0.0         0.0            0.0   0.0994115
    SiteB-4          │   0.056811         0.0            0.0   0.0571087
    SiteB-5          │        0.0   0.0128716  …    0.196689         0.0

Three methods will be demonstrated.

### Jaccard Similarity

### Czekanowski Index

First, let’s compose a syntopic table object from the “y” sample data
and extract the syntopic tables in matrix format.

``` julia
syn_y = VegSci.compose_syntopic_table_object("Sample", y)
syn_y_mat = extract_syntopic_matrix(syn_y)
syn_1_mat = extract_syntopic_matrix(syn_1)
syn_2_mat = extract_syntopic_matrix(syn_2)
```

    1×18 Named Matrix{Float64}
    A ╲ B │   Species1    Species2  …   Species19   Species20
    ──────┼──────────────────────────────────────────────────
    Syn2  │  0.0863782    0.165752  …  0.00376953    0.238411

Now we have three matrices, containg the relative frequencies of each
species present in the sample releves which constitute the
phytocoenosis’. However, each of the phytocoenosis is composed of a
different set of species, in Julia we need a helper function to merge
these matrices and ensure each matric contains each species across all
the matrices.

``` julia
merged_syn_mats = VegSci.merge_namedarrays([syn_y_mat, syn_1_mat, syn_2_mat])
```

    3×27 Named Matrix{Float64}
     A ╲ B │   Species1    Species2  …   Species15   Species18
    ───────┼──────────────────────────────────────────────────
    Sample │   0.056811   0.0128716  …         0.0         0.0
    Syn1   │   0.104218   0.0393323       0.130846   0.0828081
    Syn2   │  0.0863782    0.165752  …   0.0743653   0.0928024

``` julia
VegSci.czekanowski_index(merged_syn_mats[[:"Sample"],:], merged_syn_mats[Not(:"Sample"), :])
```

    1×2 Named Matrix{Float64}
     A ╲ B │     Syn1      Syn2
    ───────┼───────────────────
    Sample │ 0.420312  0.524003

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
