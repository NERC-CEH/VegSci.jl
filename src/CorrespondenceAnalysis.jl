using NamedArrays
using LinearAlgebra
using LinearMaps
using Arpack
using SparseArrays
using Statistics
using Random

"""
correspondence_analysis(N::NamedMatrix)::NamedMatrix

Perform a correspondence analysis following the computational algorithm outlined in appendix A of Greenacre (2017).

### Input

- `N` -- A releve by species matrix of the class NamedArrays::NamedMatrix
- `axes_n` -- The number of axes, passed to the nsv argument (Number of singular values) of `Arpack.svds`.

### Output

### Notes

### Algorithm

### References

Greenacre, M., 2017. Correspondence Analysis in Practice, Third Edition. CRC Press.

"""
function correspondence_analysis(N::NamedMatrix; axes_n::Int64 = 3)

  # N = VegSci.generate_test_array(rown = 1000, coln = 100000, meancoloccs = 15, sparse_array = true)
  # N = VegSci.generate_test_array(rown = 1000, coln = 100000, meancoloccs = 15, val_type = "bool", sparse_array = true)
  # axes_n = 3

  # Perform checks on input matrix
  if any(x -> x .< 0.0, N) == true
    println("Matrix cannot contain negative values.")
    return
  end

  if all(x -> x .== 0.0, N) == true
    println("Matrix cannot contain all zero values")
    return
  end

  # Create axes names
  axes_names = vec([string("CA")].*string.([1:1:axes_n;]))

  # A.1 Create the correspondence matrix
  P = N / sum(N)

  # A.2 Calculate column and row masses
  r = vec(sum(P, dims = 2))
  c = vec(sum(P, dims = 1))

  # A.3 Diagonal matrices of row and column masses
  Dr = sparse(NamedArray(Diagonal(r), names = (names(N)[1], names(N)[1]), dimnames = ("Rows", "Rows")))
  Dc = sparse(NamedArray(Diagonal(c), names = (names(N)[2], names(N)[2]), dimnames = ("Cols", "Cols")))

  # A.4 Calculate the matrix of standardized residuals
  SR = Dr^(-1/2) * (P - sparse(r * transpose(c))) * Dc^(-1/2)

  # A.5 Calculate the Singular Value Decomposition
  Z = Arpack.svds(LinearMap(SR), nsv = axes_n)[1]

  # A.6 Standard coordinates Φ of rows
  Φ = NamedArray(Dr^(-1/2) * Z.U, names = (names(N)[1], axes_names), dimnames = ("Rows", "Axis"))

  # A.7 Standard coordinates Γ of columns
  Γ = NamedArray(Dc^(-1/2) * Z.V, names = (names(N)[2], axes_names), dimnames = ("Cols", "Axis"))

  # A.8 Principal coordinates F of rows
  F = Φ * Diagonal(Z.S)

  # A.9 Principal coordinates G of columns
  G = Γ * Diagonal(Z.S)

  results = Dict("StandardCoordsRows" => Φ,
                 "StandardCoordsCols" => Γ,
                 "PrincipleCoordsRows" => F,
                 "PrincipleCoordsCols" => G)

  return results
  
end