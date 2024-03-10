using NamedArrays
using LinearAlgebra

"""
correspondence_analysis(N::NamedMatrix)::NamedMatrix

Perform a correspondence analysis following the computational algorithm outlined in appendix A of Greenacre (2017).

### Input

- `N` -- A releve by species matrix of the class NamedArrays::NamedMatrix

### Output

### Notes

### Algorithm

### References

Greenacre, M., 2017. Correspondence Analysis in Practice, Third Edition. CRC Press.

"""
function correspondence_analysis(N::NamedMatrix)

    # Perform checks on input matrix
    if any(x -> x .< 0.0, N) == true
      println("Matrix cannot contain negative values.")
      return
    end

    if all(x -> x .== 0.0, N) == true
      println("Matrix cannot contain all zero values")
      return
    end

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
    Φ_colnames = vec(["Dim"].*string.([1:1:size(D,1);]))
    Φ = NamedArray(Dr^(-1/2) * U, names = (Φ_rownames, Φ_colnames), dimnames = ("Plot", "Dimension"))[1:end,1:end .!= end]
    
    # A.7 Standard coordinates Γ of columns
    Γ_rownames = names(N)[2]
    Γ_colnames = vec(["Dim"].*string.([1:1:size(D,1);]))
    Γ = NamedArray(Dc^(-1/2) * V, names = (Γ_rownames, Γ_colnames), dimnames = ("Species", "Dimension"))[1:end,1:end .!= end]
    
    # A.8 Principal coordinates F of rows
    # F = Φ * D
    F = Dr^(-1/2) * U * D
    F = F[1:end,1:end .!= end]
    
    # A.9 Principal coordinates G of columns
    # G = Γ * D
    G = Dc^(-1/2) * V * D
    G = G[1:end,1:end .!= end]
  
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