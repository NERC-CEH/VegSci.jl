# Correspondence Analysis


## Outline

This page details how to perform correspondence analysis in Julia
following the computational algorithm outlined in appendix A of
Greenacre (2017) and implemented in the R package `ca` (Nenadic and
Greenacre 2007).

The Julia package MultivariateStats.jl does not currently contain an
implementation of correspondence analysis.

## Import Required Packages

## Julia

``` julia
using Pkg; Pkg.activate("docs")
using EcoVeg
using NamedArrays
using LinearAlgebra
using CSV
using BenchmarkTools
using DataFrames
# using Plots
```

## R

``` r
library(ca)
library(microbenchmark)
library(ggplot2)
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
library(JuliaCall)
```

## Create example data

Create example data in the form of a site by species named matrix, N.

<details class="code-fold">
<summary>Code</summary>

``` julia
N = generate_test_array(rown = 20, coln = 30, meancoloccs = 10, rowprefix = "SiteA-", colprefix = "Species")
```

</details>

    20×30 Named Matrix{Float64}
    Releve ╲ Species │    Species1     Species2  …    Species29    Species30
    ─────────────────┼──────────────────────────────────────────────────────
    SiteA-1          │         0.0          0.0  …          0.0          0.0
    SiteA-2          │         0.0          0.0             0.0          0.0
    SiteA-3          │   0.0240448          0.0       0.0437783          0.0
    SiteA-4          │         0.0   0.00663294             0.0    0.0297897
    SiteA-5          │         0.0          0.0       0.0390977          0.0
    SiteA-6          │         0.0          0.0       0.0742126    0.0663092
    SiteA-7          │     0.10077          0.0             0.0          0.0
    SiteA-8          │         0.0     0.100876             0.0          0.0
    ⋮                            ⋮            ⋮  ⋱            ⋮            ⋮
    SiteA-13         │    0.028682    0.0108176        0.106406          0.0
    SiteA-14         │  0.00136337     0.102277       0.0667364    0.0599434
    SiteA-15         │         0.0     0.181269             0.0     0.194195
    SiteA-16         │         0.0          0.0             0.0          0.0
    SiteA-17         │         0.0          0.0             0.0          0.0
    SiteA-18         │   0.0777982          0.0        0.103745          0.0
    SiteA-19         │         0.0          0.0             0.0          0.0
    SiteA-20         │         0.0    0.0883831  …          0.0    0.0375208

## A.1 Create the correspondence matrix

Calculate the correspondence matrix *P* following **?@eq-p**.

``` math
P = \frac{1}{n}N
```

``` julia
begin
  n = sum(N)
  P = N / n
end
```

    20×30 Named Matrix{Float64}
    Releve ╲ Species │    Species1     Species2  …    Species29    Species30
    ─────────────────┼──────────────────────────────────────────────────────
    SiteA-1          │         0.0          0.0  …          0.0          0.0
    SiteA-2          │         0.0          0.0             0.0          0.0
    SiteA-3          │  0.00120224          0.0      0.00218892          0.0
    SiteA-4          │         0.0  0.000331647             0.0   0.00148949
    SiteA-5          │         0.0          0.0      0.00195489          0.0
    SiteA-6          │         0.0          0.0      0.00371063   0.00331546
    SiteA-7          │  0.00503848          0.0             0.0          0.0
    SiteA-8          │         0.0    0.0050438             0.0          0.0
    ⋮                            ⋮            ⋮  ⋱            ⋮            ⋮
    SiteA-13         │   0.0014341   0.00054088      0.00532032          0.0
    SiteA-14         │  6.81687e-5   0.00511383      0.00333682   0.00299717
    SiteA-15         │         0.0   0.00906347             0.0   0.00970975
    SiteA-16         │         0.0          0.0             0.0          0.0
    SiteA-17         │         0.0          0.0             0.0          0.0
    SiteA-18         │  0.00388991          0.0      0.00518724          0.0
    SiteA-19         │         0.0          0.0             0.0          0.0
    SiteA-20         │         0.0   0.00441915  …          0.0   0.00187604

## A.2 Calculate column and row masses

Calculate the row and and column masses using **?@eq-row_masses** and
**?@eq-column_masses** respectively.

``` math
r = P1 \space \space
r_{i} = \sum^{J}_{j = 1} P_{ij}
```

``` math
c = P^{t}1 \space \space
c_{j} = \sum^{I}_{i = 1} P_{ij}
```

``` julia
r = vec(sum(P, dims = 2))
```

    20-element Vector{Float64}:
     0.05
     0.05
     0.05000000000000001
     0.05000000000000001
     0.05
     0.05
     0.05000000000000002
     0.049999999999999996
     0.05
     0.049999999999999996
     0.05000000000000001
     0.05
     0.050000000000000024
     0.05000000000000001
     0.05000000000000002
     0.05
     0.05
     0.05
     0.05000000000000002
     0.05000000000000001

``` julia
c = vec(sum(P, dims = 1))
```

    30-element Vector{Float64}:
     0.013569506046829185
     0.033647791130904985
     0.03988483239278359
     0.032696773065988935
     0.050093646865967255
     0.03843094071510966
     0.04510911432202968
     0.009737735123891646
     0.0534445134231606
     0.018784137551736375
     ⋮
     0.05236698046324569
     0.02790198764683744
     0.044944112895133476
     0.03048539841706622
     0.035968969779627255
     0.02513934278052018
     0.02652174455087298
     0.030160229836867728
     0.021158699754477785

## A.3 Diagonal matrices of row and column masses

``` julia
Dr = Diagonal(r)
```

    20×20 Diagonal{Float64, Vector{Float64}}:
     0.05   ⋅     ⋅     ⋅     ⋅     ⋅    …   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅    0.05   ⋅     ⋅     ⋅     ⋅        ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅    0.05   ⋅     ⋅     ⋅        ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅    0.05   ⋅     ⋅        ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅    0.05   ⋅        ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅    0.05  …   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅        ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅        ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅        ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅        ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    …   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅        ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅        ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅        ⋅     ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅       0.05   ⋅     ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    …   ⋅    0.05   ⋅     ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅        ⋅     ⋅    0.05   ⋅     ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅        ⋅     ⋅     ⋅    0.05   ⋅     ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅        ⋅     ⋅     ⋅     ⋅    0.05   ⋅ 
      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅        ⋅     ⋅     ⋅     ⋅     ⋅    0.05

``` julia
Dc = Diagonal(c)
```

    30×30 Diagonal{Float64, Vector{Float64}}:
     0.0135695   ⋅          ⋅         …   ⋅          ⋅          ⋅ 
      ⋅         0.0336478   ⋅             ⋅          ⋅          ⋅ 
      ⋅          ⋅         0.0398848      ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅             ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅             ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅         …   ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅             ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅             ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅             ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅             ⋅          ⋅          ⋅ 
     ⋮                                ⋱                        
      ⋅          ⋅          ⋅             ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅             ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅             ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅             ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅         …   ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅             ⋅          ⋅          ⋅ 
      ⋅          ⋅          ⋅            0.0265217   ⋅          ⋅ 
      ⋅          ⋅          ⋅             ⋅         0.0301602   ⋅ 
      ⋅          ⋅          ⋅             ⋅          ⋅         0.0211587

## A.4 Calculate the matrix of standardized residuals

``` math
SR = D_{r}^{-\frac{1}{2}}(P - rc^{T})D_{c}^{-\frac{1}{2}}
```

``` julia
SR = Dr^(-1/2) * (P - r * transpose(c)) * Dc^(-1/2)
```

    20×30 Named Matrix{Float64}
    Releve ╲ Species │    Species1     Species2  …    Species29    Species30
    ─────────────────┼──────────────────────────────────────────────────────
    SiteA-1          │  -0.0260476   -0.0410169  …   -0.0388331   -0.0325259
    SiteA-2          │  -0.0260476   -0.0410169      -0.0388331   -0.0325259
    SiteA-3          │    0.020108   -0.0410169       0.0175341   -0.0325259
    SiteA-4          │  -0.0260476   -0.0329313      -0.0388331    0.0132679
    SiteA-5          │  -0.0260476   -0.0410169       0.0115076   -0.0325259
    SiteA-6          │  -0.0260476   -0.0410169       0.0567201    0.0694069
    SiteA-7          │    0.167386   -0.0410169      -0.0388331   -0.0325259
    SiteA-8          │  -0.0260476    0.0819518      -0.0388331   -0.0325259
    ⋮                            ⋮            ⋮  ⋱            ⋮            ⋮
    SiteA-13         │   0.0290094   -0.0278302       0.0981715   -0.0325259
    SiteA-14         │  -0.0234305     0.083659        0.047094    0.0596213
    SiteA-15         │  -0.0260476     0.179952      -0.0388331     0.265998
    SiteA-16         │  -0.0260476   -0.0410169      -0.0388331   -0.0325259
    SiteA-17         │  -0.0260476   -0.0410169      -0.0388331   -0.0325259
    SiteA-18         │    0.123291   -0.0410169       0.0947447   -0.0325259
    SiteA-19         │  -0.0260476   -0.0410169      -0.0388331   -0.0325259
    SiteA-20         │  -0.0260476    0.0667228  …   -0.0388331    0.0251524

## A.5 Calculate the Singular Value Decomposition (SVD) of S

``` julia
begin
  svd = LinearAlgebra.svd(SR)
  U = svd.U
  S = svd.S
  V = svd.V
  Vt = svd.Vt
end
```

    20×30 Matrix{Float64}:
     -0.0148318    0.301758    0.0368141    …   0.171262    0.0305251   0.265716
      0.115241     0.0242448   0.109325        -0.117436    0.119219   -0.0164511
     -0.012284     0.131033   -0.0573014       -0.208654   -0.01347     0.168166
      0.102854    -0.191181    0.0452611        0.152765    0.14083    -0.26687
      0.0094902    0.120442    0.359512         0.170172    0.0125139   0.106037
      0.0758245   -0.0808099   0.0127464    …   0.124595   -0.0270998  -0.355176
      0.110126    -0.0410802  -0.342928        -0.200906    0.115739   -0.0707024
      0.0317845    0.206821    0.180339        -0.0393066  -0.096138    0.0299901
     -0.274413     0.178246    0.100619        -0.0660056   0.175373    0.0789571
      0.370762    -0.0926431   0.11307         -0.411728    0.136504   -0.0332505
      0.0719231    0.0206327   0.144513     …  -0.14657    -0.333402    0.11481
      0.38569     -0.122172   -0.0661283        0.233834    0.259105    0.216812
      0.0162697   -0.0862636   0.041274         0.339983   -0.157945   -0.154871
      0.142567     0.12103    -0.201633        -0.0173558   0.174748    0.0407328
      0.105553    -0.319151   -0.00966001      -0.0193188   0.28795     0.177635
      0.00120994  -0.0995065  -0.00597018   …  -0.295173   -0.149831    0.0797212
     -0.129409    -0.288619    0.000480026      0.158042   -0.0673486   0.0467113
     -0.241891    -0.398829    0.0848882       -0.195579    0.101503    0.130363
     -0.00252854   0.015      -0.148526        -0.0979839  -0.0764976   0.285176
      0.198088     0.0180075   0.328286         0.0976223   0.0140308   0.243499

## A.6 Standard coordinates *Φ* of rows

``` math
\Phi = D_{r}^{-\frac{1}{2}} U
```

``` julia
Φ = Dr^(-1/2) * U
```

    20×20 Matrix{Float64}:
      0.242949  -0.894712    1.13509    …   0.00259932   0.497571    -1.0
     -1.92928   -1.0432      0.264968       0.390087    -0.587302    -1.0
      0.521654  -0.167549   -1.98314       -0.544606    -0.637304    -1.0
     -1.83683   -1.2923      0.842607       0.22697     -1.00396     -1.0
      1.08557   -0.789531   -2.06548        1.25621     -0.358405    -1.0
      0.422749   0.408169    0.815066   …   1.35711     -1.41815     -1.0
      0.291823  -0.350631    0.658339      -0.992634    -0.541535    -1.0
      0.202988   1.1853     -0.256498      -1.56327      0.482102    -1.0
      0.312707  -0.651128   -0.029152       1.96802      0.466441    -1.0
     -0.340455  -0.380844    0.163086      -0.635481     3.0073      -1.0
     -0.462721  -0.57005    -1.14467    …  -0.0993965   -0.0136879   -1.0
      1.22569    0.323578    1.88322        0.109187    -0.797712    -1.0
     -0.304117   0.0566031   0.538648      -0.723149     0.246609    -1.0
      0.307342   0.108982   -0.0143044     -1.86646     -1.25533     -1.0
      2.35511    0.066256    1.07961        0.731087     0.927145    -1.0
     -1.10656    2.21582     0.235054   …   0.0124529    0.75946     -1.0
     -0.55389   -0.263702   -0.839243       0.820769     1.31725     -1.0
     -0.53447    2.77323    -0.714349       1.07913     -0.67428     -1.0
     -0.845596  -0.105057    0.346835      -0.130168    -0.00721748  -1.0
      0.94535   -0.629225   -0.91568       -1.39845     -0.408986    -1.0

``` julia
# NamedArrays.setnames!(NamedArray(Φ), names = vec(names(N)[1]))
# NamedArray(Φ, names(N)[1])
```

## A.7 Standard coordinates *Γ* of columns

``` math
\Gamma = D_{c}^{-\frac{1}{2}} V
```

``` julia
Γ = Dc^(-1/2) * V
```

    30×20 Matrix{Float64}:
     -0.127325   0.989289   -0.105453   …  -2.07653    -0.0217064   1.7005
      1.64506    0.132172    0.714333      -2.17424     0.0817734   0.098169
      0.184336   0.547411   -0.28692        0.425053   -0.7437      1.6438
     -0.42474    1.67927    -1.13484       -0.292643    0.664241    0.630497
      1.57816    0.247273    1.61133        0.570176   -0.476515    0.360303
     -0.990649  -0.433224    0.101574   …  -2.05811    -0.0110217   0.531588
      0.400289  -0.398498   -1.10858       -1.01338     1.32371     1.21158
     -0.199446  -0.0224278  -0.459178       2.78177     0.840407   -0.2313
     -0.752107  -1.00735    -0.47954        0.457041   -0.374511    0.832762
      0.957237  -0.705942   -2.59576        2.34727    -0.568138    0.711355
      ⋮                                 ⋱                          
     -0.647238  -0.570337   -0.724105      -0.663303    0.999029    1.13825
     -0.150498  -0.508282    1.23055        1.6668      1.28966     0.784872
      0.329692  -0.710874   -0.110099      -0.0776507  -0.746994    0.856759
     -0.918172  -0.780123    0.557805      -0.819415    1.19215     0.907027
     -0.45945   -0.753904   -1.52835    …   0.532503   -0.585457    1.1874
      0.192537   0.601972    0.321879      -1.14566    -1.77827     1.47708
      1.05162   -0.72111    -1.28122       -1.20094    -0.601664    0.599443
      0.175768   0.686478   -0.0775625      0.584472   -0.440485    0.0807912
      1.82673   -0.113097    1.1561         0.896212    1.96051     1.67399

## A.8 Principal coordinates F of rows

``` math
F = D_{r}^{-\frac{1}{2}} U D_{\alpha} = \Phi D_{\alpha}
```

``` julia
F = Φ * Diagonal(S)
```

    20×20 Matrix{Float64}:
      0.148499  -0.544564    0.601765    …   0.0431903    -1.50475e-16
     -1.17925   -0.634941    0.140471       -0.0509791    -1.50475e-16
      0.318854  -0.101978   -1.05135        -0.0553194    -1.50475e-16
     -1.12273   -0.786556    0.446704       -0.0871461    -1.50475e-16
      0.663536  -0.480546   -1.09501        -0.0311104    -1.50475e-16
      0.258399   0.248431    0.432103    …  -0.123099     -1.50475e-16
      0.178372  -0.21341     0.349015       -0.0470064    -1.50475e-16
      0.124074   0.721429   -0.135981        0.0418475    -1.50475e-16
      0.191138  -0.396307   -0.0154548       0.0404881    -1.50475e-16
     -0.208098  -0.2318      0.0864593       0.26104      -1.50475e-16
     -0.282832  -0.346959   -0.606839    …  -0.00118814   -1.50475e-16
      0.749184   0.196945    0.998378       -0.0692431    -1.50475e-16
     -0.185887   0.0344513   0.285562        0.0214062    -1.50475e-16
      0.187859   0.0663314  -0.00758342     -0.108965     -1.50475e-16
      1.43953    0.0403265   0.572348        0.0804782    -1.50475e-16
     -0.67637    1.34865     0.124613    …   0.0659227    -1.50475e-16
     -0.338557  -0.160502   -0.44492         0.11434      -1.50475e-16
     -0.326687   1.68792    -0.378709       -0.058529     -1.50475e-16
     -0.516859  -0.0639427   0.183873       -0.000626493  -1.50475e-16
      0.577832  -0.382976   -0.485443       -0.0355009    -1.50475e-16

## A.9 Principal coordinates G of columns

``` math
G = D_{c}^{-\frac{1}{2}} V D_{\alpha} = \Gamma D_{\alpha}
```

``` julia
G = Γ * Diagonal(S)
```

    30×20 Matrix{Float64}:
     -0.0778254   0.602128   -0.0559054  …  -0.00188417    2.55882e-16
      1.00552     0.0804463   0.3787         0.00709811    1.4772e-17
      0.112673    0.33318    -0.152109      -0.0645547     2.4735e-16
     -0.259616    1.02209    -0.601632       0.0576576     9.48738e-17
      0.96463     0.150502    0.85424       -0.0413625     5.42165e-17
     -0.60552    -0.263681    0.0538489  …  -0.000956704   7.99906e-17
      0.244671   -0.242545   -0.587711       0.114901      1.82311e-16
     -0.121908   -0.0136506  -0.243431       0.0729492    -3.48049e-17
     -0.459715   -0.613121   -0.254226      -0.0325084     1.2531e-16
      0.585097   -0.42967    -1.37613       -0.0493156     1.07041e-16
      ⋮                                  ⋱                
     -0.395615   -0.347134   -0.383881       0.0867179     1.71278e-16
     -0.0919896  -0.309364    0.652373       0.111946      1.18103e-16
      0.201519   -0.432672   -0.0583683     -0.0648407     1.28921e-16
     -0.56122    -0.47482     0.295718       0.103481      1.36485e-16
     -0.280832   -0.458862   -0.810247   …  -0.050819      1.78674e-16
      0.117686    0.366389    0.170643      -0.154357      2.22263e-16
      0.642789   -0.438902   -0.679235      -0.0522257     9.0201e-17
      0.107436    0.417823   -0.0411194     -0.038235      1.2157e-17
      1.11656    -0.068836    0.612899       0.170177      2.51893e-16

## A.10 Principal inertias *λ*<sub>*k*</sub>

``` math
\lambda_{k} = \alpha_{k}^{2}, k = 1,2,...,\space where \space k = min\{I-1,J-1\}
```

``` julia
F * Dr * transpose(F)
```

    20×20 Matrix{Float64}:
      0.133386     -0.0199282   -0.00262721  …  -0.000930025   0.00522482
     -0.0199282     0.165131    -0.00387565      1.64847e-6   -0.021691
     -0.00262721   -0.00387565   0.163622       -0.0180622     0.0014463
      0.0136202     0.0724153   -0.0361814       0.0481068    -0.0373496
     -0.00653076   -0.05         0.0492494      -0.0274908     0.0375868
      0.000510675  -0.0187512   -0.0275665   …  -0.0329749    -0.00489999
      0.0211456     0.00469376  -0.0223489      -0.0277637    -0.0176429
     -0.0394234    -0.0203218   -0.00893078     -0.0098347    -0.0197946
      0.00955613    0.00749034   0.00653559      0.00928847    0.0361745
     -0.0130638     0.0209051    0.0114159       0.0125581    -0.0112933
     -0.0372901     0.0129624   -0.0428679   …  -0.00196722   -0.0117365
      0.0256133    -0.0447863   -0.0357827      -0.0129672    -0.0153736
      0.0378716    -0.0146924   -0.0155792      -0.00845823   -0.00334006
     -0.0241973    -0.0148991    0.00101534     -0.0195392    -0.030539
     -0.0200104    -0.05        -0.00722887     -0.0223349     0.00815378
     -0.0258622     0.0108922   -0.0307117   …  -0.0176206    -0.0369671
     -0.00833557    0.0144544    0.00727858     -0.0386094    -0.00959083
     -0.048729     -0.05         0.0111995       0.0271367    -0.0199467
     -0.000930025   1.64847e-6  -0.0180622       0.140567      0.000894608
      0.00522482   -0.021691     0.0014463       0.000894608   0.150684

``` julia
G * Dr * transpose(G)
```

    30×30 Matrix{Float64}:
      0.20946     -0.0475376   -0.0019206   …   0.0463967   -0.0373442
     -0.0475376    0.117078     0.0264673      -0.0140239    0.107478
     -0.0019206    0.0264673    0.0900501      -0.00920603   0.0179431
      0.0854656   -0.0304974    0.0624372       0.0430568   -0.0399154
      0.00996817   0.0619846   -0.00426866     -0.0194572    0.0951109
     -0.0354523   -0.0213607    0.0153405   …  -0.0159621   -0.00963823
     -0.0214092   -0.0277921   -0.0224559      -0.00834068  -0.036177
     -0.0489229   -0.0174135   -0.0397992       0.00487785   0.0200531
      0.0125807   -0.0178447   -0.00709734     -0.0347688   -0.038391
     -0.00493186  -0.0190991   -0.0196862       0.0166848   -0.0388197
      ⋮                                     ⋱               
     -0.0213967   -0.00326276  -0.0296029      -0.00142531  -0.0184714
     -0.0301548   -0.0100509    0.0345669       0.00257554  -0.0483139
      0.0200685   -0.02622     -0.0227733      -0.00921205  -0.0108574
      0.032985    -0.0256453   -0.00478235     -0.00765175  -0.0156708
     -0.0495778   -0.0143627   -0.010305    …  -0.02937     -0.0282779
     -0.0218243    0.00860405  -0.00304252      0.0403398   -0.000514609
     -0.0320589    0.00894107   0.0298837      -0.00363707  -0.014047
      0.0463967   -0.0140239   -0.00920603      0.0788655   -0.00246307
     -0.0373442    0.107478     0.0179431      -0.00246307   0.22503

## Create Correspondence Analysis Function

``` julia
function correspondence_analysis(N::NamedMatrix)
  
  # A.1 Create the correspondence matrix
  P = N / sum(N)

  # A.2 Calculate column and row masses
  r = vec(sum(P, dims = 2))
  c = vec(sum(P, dims = 1))

  # A.3 Diagonal matrices of row and column masses
  Dr = Diagonal(r)
  Dc = Diagonal(c)

  # A.4 Calculate the matrix of standardized residuals
  SR = Dr^(-1/2) * (P - r * transpose(c)) * Dc^(-1/2)

  # A.5 Calculate the Singular Value Decomposition (SVD) of S
  svd = LinearAlgebra.svd(SR)
  U = svd.U
  V = svd.V
  S = svd.S
  D = Diagonal(S)

  # A.6 Standard coordinates Φ of rows
  Φ_rownames = names(N)[1]
  Φ_colnames = vec(["Dim"].*string.([1:1:size(N,1);]))
  Φ = NamedArray(Dr^(-1/2) * U, names = (Φ_rownames, Φ_colnames), dimnames = ("Plot", "Dimension"))
  
  # A.7 Standard coordinates Γ of columns
  Γ_rownames = names(N)[2]
  Γ_colnames = vec(["Dim"].*string.([1:1:size(N,1);]))
  Γ = NamedArray(Dc^(-1/2) * V, names = (Γ_rownames, Γ_colnames), dimnames = ("Species", "Dimension"))
  
  # A.8 Principal coordinates F of rows
  # F = Φ * D
  F = Dr^(-1/2) * U * D
  
  # A.9 Principal coordinates G of columns
  # G = Γ * D
  G = Dc^(-1/2) * V * D

  # [1:end, 1:end .∉ [[20]]]

  results = (sv = D, # Singular values
             rownames = names(N)[1], # Row names
             rowmass = r, # Row masses
            #  rowdist = , # Row chi-square distances to centroid
            #  rowinertia = , # Row inertias
             rowcoord = Φ, # Row standard coordinates
            #  rowsup = , # Indicies of row supplementary points
             colnames = names(N)[2], # Column names
             colmass = c, # Column masses
            #  coldist = , # Column chi-square distances to centroid
            #  colinertia = , # Column inertias
             colcoord = Γ, # Column standard coordinates
            #  colsup = , # Indices of column supplementary points
            N = N # The frequency table
            )

  return results

end
```

    correspondence_analysis (generic function with 1 method)

### Test Function

``` julia
ca_results = correspondence_analysis(N)
```

    (sv = [0.6112355805256184 0.0 … 0.0 0.0; 0.0 0.6086472613581309 … 0.0 0.0; … ; 0.0 0.0 … 0.0868021663091435 0.0; 0.0 0.0 … 0.0 1.5047469066961423e-16], rownames = ["SiteA-1", "SiteA-2", "SiteA-3", "SiteA-4", "SiteA-5", "SiteA-6", "SiteA-7", "SiteA-8", "SiteA-9", "SiteA-10", "SiteA-11", "SiteA-12", "SiteA-13", "SiteA-14", "SiteA-15", "SiteA-16", "SiteA-17", "SiteA-18", "SiteA-19", "SiteA-20"], rowmass = [0.05, 0.05, 0.05000000000000001, 0.05000000000000001, 0.05, 0.05, 0.05000000000000002, 0.049999999999999996, 0.05, 0.049999999999999996, 0.05000000000000001, 0.05, 0.050000000000000024, 0.05000000000000001, 0.05000000000000002, 0.05, 0.05, 0.05, 0.05000000000000002, 0.05000000000000001], rowcoord = [0.24294884884619708 -0.8947124099801201 … 0.4975713851403658 -0.9999999999999994; -1.9292845375530772 -1.0431998727739231 … -0.5873018014013842 -1.0000000000000004; … ; -0.8455963688078192 -0.10505701334704681 … -0.00721748166458637 -0.9999999999999997; 0.9453499670603035 -0.6292248358098728 … -0.40898597768133454 -1.0000000000000007], colnames = ["Species1", "Species2", "Species3", "Species4", "Species5", "Species6", "Species7", "Species8", "Species9", "Species10"  …  "Species21", "Species22", "Species23", "Species24", "Species25", "Species26", "Species27", "Species28", "Species29", "Species30"], colmass = [0.013569506046829185, 0.033647791130904985, 0.03988483239278359, 0.032696773065988935, 0.050093646865967255, 0.03843094071510966, 0.04510911432202968, 0.009737735123891646, 0.0534445134231606, 0.018784137551736375  …  0.04415719303484158, 0.05236698046324569, 0.02790198764683744, 0.044944112895133476, 0.03048539841706622, 0.035968969779627255, 0.02513934278052018, 0.02652174455087298, 0.030160229836867728, 0.021158699754477785], colcoord = [-0.12732473518133222 0.9892893719773834 … -0.021706430526786866 1.7004961176270454; 1.645056250558114 0.13217224980588455 … 0.08177343835660492 0.09816901253790525; … ; 0.17576777088639256 0.6864779779908547 … -0.4404845840447725 0.08079116906761726; 1.8267280556600065 -0.11309667653297997 … 1.9605105724362382 1.6739874913739352], N = [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.08838309173227708 … 0.0 0.0375208103875558])

### Plot Test Results

``` julia
begin
  plots_x = ca_results.rowcoord[:,"Dim1"]
  plots_y = ca_results.rowcoord[:,"Dim2"]

  species_x = ca_results.colcoord[:,"Dim1"]
  species_y = ca_results.colcoord[:,"Dim2"]

  # scatter(species_x, species_y, series_annotations = text.(ca_results.colnames, 8, :right, :bottom), framestyle=:origin, label = "Species")
  # scatter!(plots_x, plots_y, series_annotations = text.(ca_results.rownames, 8, :right, :bottom), label = "Plots")
end
```

    30-element Named Vector{Float64}
    Species   │ 
    ──────────┼───────────
    Species1  │   0.989289
    Species2  │   0.132172
    Species3  │   0.547411
    Species4  │    1.67927
    Species5  │   0.247273
    Species6  │  -0.433224
    Species7  │  -0.398498
    Species8  │ -0.0224278
    ⋮                    ⋮
    Species23 │  -0.508282
    Species24 │  -0.710874
    Species25 │  -0.780123
    Species26 │  -0.753904
    Species27 │   0.601972
    Species28 │   -0.72111
    Species29 │   0.686478
    Species30 │  -0.113097

## R vs Julia Dune Data

Let’s see how the R function `ca::ca()` compares to the Julia function
`correspondence_analysis()` defined above.

First load the `dune` dataset bundled in `{vegan}`.

## Julia

``` julia
begin
    dune_df = CSV.read("../data/dune.csv", DataFrame, header = 1)
    dune_na = NamedArray(Array(dune_df))
    NamedArrays.setnames!(dune_na, names(dune_df), 2)
end
```

    (OrderedCollections.OrderedDict{Any, Int64}("1" => 1, "2" => 2, "3" => 3, "4" => 4, "5" => 5, "6" => 6, "7" => 7, "8" => 8, "9" => 9, "10" => 10…), OrderedCollections.OrderedDict{Any, Int64}("Achimill" => 1, "Agrostol" => 2, "Airaprae" => 3, "Alopgeni" => 4, "Anthodor" => 5, "Bellpere" => 6, "Bromhord" => 7, "Chenalbu" => 8, "Cirsarve" => 9, "Comapalu" => 10…))

## R

``` r
dune <- as.matrix(read.csv(file = "../data/dune.csv"))
```

### Run CA in Julia on Dune data

``` julia
dune_ca_julia = correspondence_analysis(dune_na)
```

    (sv = [0.7321237072731497 0.0 … 0.0 0.0; 0.0 0.6325690627910964 … 0.0 0.0; … ; 0.0 0.0 … 0.05896537573383853 0.0; 0.0 0.0 … 0.0 7.091790502929597e-17], rownames = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"], rowmass = [0.02627737226277372, 0.06131386861313868, 0.05839416058394161, 0.0656934306569343, 0.06277372262773723, 0.07007299270072993, 0.05839416058394161, 0.05839416058394161, 0.06131386861313869, 0.06277372262773723, 0.04671532846715329, 0.051094890510948905, 0.04817518248175183, 0.03503649635036496, 0.03357664233576642, 0.04817518248175183, 0.0218978102189781, 0.03941605839416058, 0.04525547445255475, 0.04525547445255475], rowcoord = [0.8116737226994164 1.0826713631197626 … 2.417299499470586 0.9999999999999976; 0.6326772260449806 0.6958357017624843 … -0.049159624022002205 0.9999999999999993; … ; 0.6902687683328907 -3.2642026272615623 … 0.014843350243563693 1.0; -1.9443804612344018 -1.0688808878703615 … -1.087127732821526 1.0000000000000004], colnames = ["Achimill", "Agrostol", "Airaprae", "Alopgeni", "Anthodor", "Bellpere", "Bromhord", "Chenalbu", "Cirsarve", "Comapalu"  …  "Ranuflam", "Rumeacet", "Sagiproc", "Salirepe", "Scorautu", "Trifprat", "Trifrepe", "Vicilath", "Bracruta", "Callcusp"], colmass = [0.02335766423357664, 0.07007299270072993, 0.0072992700729927005, 0.05255474452554745, 0.030656934306569343, 0.01897810218978102, 0.021897810218978103, 0.00145985401459854, 0.00291970802919708, 0.00583941605839416  …  0.020437956204379562, 0.026277372262773723, 0.029197080291970802, 0.016058394160583942, 0.07883211678832117, 0.01313868613138686, 0.06861313868613139, 0.00583941605839416, 0.07153284671532847, 0.014598540145985401], colcoord = [1.2410388912758308 -0.13374974813501736 … 1.0922556626083908 0.5378610945054166; -1.275443784762892 0.3264688587581394 … -1.3054088243241777 0.6482522649395861; … ; -0.2488978439865396 -0.4185696756019523 … -0.01770675187655213 -0.49012289845882706; -2.6661966956330643 -0.8970175654454415 … -0.22298366681774884 -0.2585309105382417], N = [1 0 … 0 0; 3 0 … 0 0; … ; 0 0 … 3 0; 0 5 … 4 3])

### Run CA in R using Dune data

``` r
dune_ca_r <- ca::ca(dune)
```

``` r
dune_ca_coords_r <- ca::cacoord(dune_ca_r)
dune_ca_r_rowcoord <- dune_ca_coords_r$rows
dune_ca_r_colcoord <- dune_ca_coords_r$columns
```

### Check equality of R and Julia CA results

Check whether the standard row (plot) and column (species) coordinates
produced by the Julia function `correspondence_analaysis()` the R
function `ca::ca()` are identical.

First retrieve the Julia results in R.

``` r
dune_ca_julia <- JuliaCall::julia_eval("dune_ca_julia")
dune_ca_julia_rowcoord <- dune_ca_julia$rowcoord
dune_ca_julia_colcoord <- dune_ca_julia$colcoord
```

Then check for equivalence. Note that there are two differences in R and
Julia data structure:

1.  Row and column names are stripped from the Julia results when
    calling `JuliaCall::julia_eval("dune_ca_julia")`, so the R results
    are also unnamed.
2.  The Julia results contain an additional dimension, populated with
    1.0 values.

These differences are corrected below whilst checking for equivalence.

The R function `all.equal` is used rather than `identical` as the Julia
implementation returns results with higher precision.

``` r
plots_allequal <- all.equal(dune_ca_julia_rowcoord[,-20], unname(dune_ca_r_rowcoord))
species_allequal <- all.equal(dune_ca_julia_colcoord[,-20], unname(dune_ca_r_colcoord))
isTRUE(all(plots_allequal, species_allequal))
```

    [1] TRUE

Let’s view the first two dimensions of the Julia and R results
side-by-side.

### Visualise standard coordinates

![](CorrespondenceAnalysis.markdown_strict_files/figure-markdown_strict/plot_dune_julia_vs_r-1.png)

### Benchmark Functions

``` julia
@benchmark correspondence_analysis(dune_na)
```

    BenchmarkTools.Trial: 10000 samples with 1 evaluation.
     Range (min … max):  106.235 μs …  2.049 ms  ┊ GC (min … max): 0.00% … 84.86%
     Time  (median):     109.876 μs              ┊ GC (median):    0.00%
     Time  (mean ± σ):   115.207 μs ± 79.273 μs  ┊ GC (mean ± σ):  2.69% ±  3.73%

      ▅█▇█▇▅▄▂▁                                                    ▂
      ██████████▇▆▅▄▅▄▆▅▆▄▅▅▅▆▇▇▇▇▇▅▅▄▄▄▄▄▄▅▄▆▄▄▆▅▅▆▆▄▅▄▅▄▄▅▄▃▄▄▄▄ █
      106 μs        Histogram: log(frequency) by time       172 μs <

     Memory estimate: 88.81 KiB, allocs estimate: 316.

``` r
microbenchmark::microbenchmark(ca::ca(dune))
```

    Unit: microseconds
             expr     min       lq     mean   median       uq     max neval
     ca::ca(dune) 444.605 450.4475 458.1625 455.0155 459.2335 582.436   100

## References

Greenacre, Michael. 2017. *Correspondence Analysis in Practice, Third
Edition*. CRC Press.

Nenadic, Oleg, and Michael Greenacre. 2007. “Correspondence Analysis in
R, with Two- and <span class="nocase">Three-dimensional Graphics</span>:
The Ca Package.” *Journal of Statistical Software* 20 (February): 1–13.
<https://doi.org/10.18637/jss.v020.i03>.
