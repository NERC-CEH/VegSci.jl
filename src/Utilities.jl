using NamedArrays
using SparseArrays
using InvertedIndices
using DataFrames

"""
    generate_test_array(;rown::Int64, coln::Int64,
                        meancoloccs::Int64,
                        rowprefix::String = "Releve", colprefix::String = "Species",
                        rowdim::String = "Releve", coldim::String = "Species",
                        sparse_array::Bool = false)

Create a `rown` by `coln` NamedArray object containing random values.

...
# Arguments
- `rown::Int64`: The number of rows in the array.
- `coln::Int64`: The number of columns in the array.
- `meancoloccs::Int64`: The mean number of non-zero elements in each column, usually representing the species richness of each Releve.
- `rowprefix::String`: The prefix to the row number. "Releve" by default.
- `colprefix::String`: The prefix to the column number. "Species" by default.
- `rowdim::String`: The row dimension name. "Releve" by default.
- `coldim::String`: The column dimension name. "Species" by default.
- `sparse_array::Bool`: If true a names sparse matrix is returned. If false a named dense matrix is returned. false by default.
...

# Examples
```julia
julia>VegSci.generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "Releve", colprefix = "Species")
```
"""
function generate_test_array(;rown::Int64, coln::Int64,
                             meancoloccs::Int64,
                             rowprefix::String = "Releve", colprefix::String = "Species",
                             rowdim::String = "Releve", coldim::String = "Species",
                             sparse_array::Bool = false)

    nonzerop = meancoloccs / coln
    rownames = vec([string("$rowprefix")].*string.([1:1:rown;]))
    colnames = vec([string("$colprefix")].*string.([1:1:coln;]))

    if sparse_array == true
        vals = sparse(sprand(Float64, rown, coln, nonzerop))
    elseif sparse_array == false
        vals = Array(sprand(Float64, rown, coln, nonzerop))
    end

    x = NamedArrays.NamedArray(vals, names = (rownames, colnames), dimnames = (rowdim, coldim))
    y = x ./ sum(x, dims = 2)

    return y

end

"""
    nm_to_df(nm::NamedMatrix)

Convert a named matrix of class NamedArrays.NamedMatrix to a data frame of class DataFrames.DataFrame

...
# Arguments
- `nm::NamedMatrix`: A named matrix of class NamedArrays.NamedMatrix
...

...
# Returns
A data frame of class DataFrames.DataFrame with column names equal to the column names 
of the named matrix and with the first column of the data frame equal to the rownames of the named matrix.
...

# Examples
```julia
julia>nm = VegSci.generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "Releve", colprefix = "Species")
julia>VegSci.nm_to_df(nm)
```
"""
function nm_to_df(nm::NamedMatrix)

    df = DataFrame(nm, Symbol.(names(nm)[2]))
    insertcols!(df, 1, Symbol(dimnames(nm, 1)) => names(nm)[1])

  return df

end

"""
Draws heavily from the function outlined here: https://discourse.julialang.org/t/nanmean-options/4994/17
"""
function nzfunc(f::Function, x::NamedArray; dims::Int64 = 1)

    _nzfunc(fn, A, ::Colon) = fn(filter(!iszero, A))
    _nzfunc(fn, A, dims) = mapslices(a -> _nzfunc(fn, a , :), A, dims = dims)
    nzfunc(fn, A; dims = :) = _nzfunc(fn, A, dims)
    y = nzfunc(f, x, dims = dims)
    setnames!(y, names(x)[2], 2)
    
    return y

end

function merge_namedarrays(mats::Vector)

    df_all = DataFrame()
    
    for mat in mats
        rowdimname = dimnames(mat, 1)
        df = DataFrame(mat, Symbol.(names(mat)[2]))
        insertcols!(df, 1, Symbol(rowdimname) => names(mat)[1])
        df = stack(df, Not(rowdimname))
        df_all = [df_all; df]
    end

    df_all_wide = unstack(df_all, fill = 0.0)
    df_all_wide_prepped = select(df_all_wide, Not([1]))

    rownames = df_all_wide[!, 1]
    colnames = names(df_all_wide_prepped)

    results = NamedArray(Matrix(df_all_wide_prepped), names = (rownames, colnames))

    return results

end

function align_array_columns(x::NamedArray, y::NamedArray, colorder::String = "x")

    # x = VegSci.generate_test_array(rown = 100, coln = 100, meancoloccs = 5, rowprefix = "SiteA-", colprefix = "Species", sparse_array = true)
    # y = VegSci.generate_test_array(rown = 50, coln = 100, meancoloccs = 5, rowprefix = "SiteB-", colprefix = "Species", sparse_array = true)

    # x = x[:,Not(["Species3", "Species10"])]
    # y = y[:,Not(["Species4"])]

    # typeof(x)
    # typeof(y)

    # Check which columns are missing from x and y
    x_missing_cols = setdiff(Set(names(y)[2]), Set(names(x)[2]))
    y_missing_cols = setdiff(Set(names(x)[2]), Set(names(y)[2]))

    x_mat = copy(x)

    # If there are missing columns in the x matrix
    if length(x_missing_cols) != 0
        x_mat_missing = NamedArray(zeros(size(x,1), length(x_missing_cols)), names = (vec(names(x)[1]), collect(x_missing_cols)))
        x_mat_colnames = names(x)[2]
        x_mat = [x x_mat_missing]
        setnames!(x_mat, [x_mat_colnames; collect(x_missing_cols)], 2)
    end

    y_mat = copy(y)

    # If there are missing columns in the x matrix
    if length(y_missing_cols) != 0
        y_mat_missing = NamedArray(zeros(size(y,1), length(y_missing_cols)), names = (vec(names(y)[1]), collect(y_missing_cols)))
        y_mat_colnames = names(y)[2]
        y_mat = [y y_mat_missing]
        setnames!(y_mat, [y_mat_colnames; collect(y_missing_cols)], 2)
    end

    if colorder == "x"
        y_mat = y_mat[:, names(x_mat)[2]]
    elseif colorder == "y"
        x_mat = x_mat[:, names(y_mat)[2]]
    end

    aligned_mats = (x = x_mat, y = y_mat)

    return aligned_mats

end