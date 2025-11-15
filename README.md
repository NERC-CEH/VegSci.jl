

![version](https://img.shields.io/badge/version-0.1.0-blue)
[![CI](https://github.com/ZekeMarshall/VegSci.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/ZekeMarshall/VegSci.jl/actions/workflows/ci.yml?query=branch%3Amain)
[![Codecov](https://codecov.io/gh/ZekeMarshall/VegSci.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ZekeMarshall/VegSci.jl)
[![Aqua
QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/ZekeMarshall/VegSci.jl)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999.png)]()

# VegSci.jl

`VegSci.jl` is a package containing tools for vegetation science using
Julia (Bezanson et al. 2017), a growing scientific programming language
which solves the 'two language problem' (Roesch et al. 2023), offering C
and FORTRAN-like performance alongside the readability and
user-friendliness of higher level languages such as Python and R. `VegSci.jl`
aims to impliment selected functionality found in popular vegetation
science software programs/packages such as JUICE (Tichý 2002), vegan
(Oksanen et al. 2022), ade4 (Dray and Dufour 2007), indicspecies (De
Cáceres and Legendre, 2009), vegclust (De
Cáceres, Font, and Oliva 2010; De Cáceres, Legendre, and He 2013), and
vegsoup (Kaiser 2021)  on an ad-hoc basis, as and when required by the package authors,
with a user-friendly API and transparent methodologies. `VegSci.jl` is being developed with the
aim of assisting in the creation of high-performance, reproducible
analytical pipelines in vegetation research (Sperandii et al. 2024).

## Design

`VegSci.jl` is structured around the Julia package `NamedArrays.jl`
facilitating the use of named relevé by species matrices, which
constitute the basic input for most functions.

## Development Notes

<details>
<summary><b>Dependency Management</b></summary>

Coming soon...

</details>

<details>
<summary><b>Version Control</b></summary>

Coming soon...

</details>

<details>
<summary><b>Release Steps</b></summary>

Coming soon...

</details>

<details>
<summary><b>Testing</b></summary>

To test `VegSci.jl`
    1.  Run `] activate .` to activate the package environment.
    2.  Run `] test`
    2.  Run `] activate` to exit the package environment.

</details>

<details>
<summary><b>Documentation</b></summary>

Run `quarto check` in the terminal to check whether Quarto is ok.

If Julia has been updated you will need to re-install the IJulia kernel
by running `using IJulia` then `installkernel("Julia")` in the Julia
terminal.

To render the README run `quarto render README.qmd --to md`

</details>

## References

Bezanson, Jeff, Alan Edelman, Stefan Karpinksi, and Viral B. Shah. 2017.
"Julia: A Fresh Approach to Numerical Computing." *SIAM Review* 59 (1):
65-98. <https://doi.org/10.1137/141000671>.

De Cáceres, Miquel, Pierre Legendre, 2009. 
Associations between species and groups of sites: 
indices and statistical inference. Ecology 90, 3566-3574. 
https://doi.org/10.1890/08-1823.1

De Cáceres, Miquel, Xavier Font, and Francesc Oliva. 2010. "The
Management of Vegetation Classifications with Fuzzy Clustering."
*Journal of Vegetation Science* 21 (6): 1138-51.
<https://doi.org/10.1111/j.1654-1103.2010.01211.x>.

De Cáceres, Miquel, Pierre Legendre, and Fangliang He. 2013.
"Dissimilarity Measurements and the Size Structure of Ecological
Communities." *Methods in Ecology and Evolution* 4 (12): 1167-77.
<https://doi.org/10.1111/2041-210X.12116>.

Dray, Stéphane, and Anne-Béatrice Dufour. 2007. "The Ade4 Package:
Implementing the Duality Diagram for Ecologists." *Journal of
Statistical Software* 22 (September): 1-20.
<https://doi.org/10.18637/jss.v022.i04>.

Kaiser, Roland. 2021. "Vegsoup: Classes and Methods for Phytosociology."

Oksanen, Jari, Gavin Simpson, Peter Solymos, Leo Lahti, Geoffrey
Hannigan, James Weedon, Eduard Szöcs, et al. 2022. "Vegandevs/Vegan:
Vegan 2.6-4 on CRAN." <https://doi.org/10.5281/zenodo.7185692>.

Roesch, Elisabeth, Joe G. Greener, Adam L. MacLean, Huda Nassar,
Christopher Rackauckas, Timothy E. Holy, and Michael P. H. Stumpf. 2023.
"Julia for Biologists." *Nature Methods* 20 (5): 655-64.
<https://doi.org/10.1038/s41592-023-01832-z>.

Sperandii, Marta Gaia, Manuele Bazzichetto, Glenda Mendieta-Leiva,
Sebastian Schmidtlein, Michael Bott, Renato A. Ferreira de Lima, Valério
D. Pillar, Jodi N. Price, Viktoria Wagner, and Milan Chytrý. 2024.
"Towards More Reproducibility in Vegetation Research." *Journal of
Vegetation Science* 35 (1): e13224. <https://doi.org/10.1111/jvs.13224>.

Tichý, Lubomír. 2002. "JUICE, Software for Vegetation Classification."
*Journal of Vegetation Science* 13 (3): 451-53.
<https://doi.org/10.1111/j.1654-1103.2002.tb02069.x>.
