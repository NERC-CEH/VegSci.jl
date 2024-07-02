using NamedArrays
using SparseArrays
using InvertedIndices
using DataFrames

"""
    nm_to_df(nm::NamedMatrix)

Convert a named matrix of class NamedArrays.NamedMatrix to a 

...
# Arguments
- `nm::NamedMatrix`: A named matric of class NamedArrays.NamedMatrix
...

...
# Returns
A data frame of class DataFrames.DataFrame with column names equal to the column names 
of the named matrix and with the first column of the data frame equal to the rownames of the named matrix.
...

# Examples
```julia
julia>generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "Releve", colprefix = "Species")
```
"""
function nm_to_df(nm::NamedMatrix)

    df = DataFrame(nm, Symbol.(names(nm)[2]))
    insertcols!(df, 1, Symbol(dimnames(nm, 1)) => names(nm)[1])

  return df

end

"""
    generate_test_array(;rown::Int64, coln::Int64,
                        meancoloccs::Int64,
                        rowprefix::String = "Releve", colprefix::String = "Species",
                        rowdim::String = "Releve", coldim::String = "Species")

Create a `rown` by `coln` NamedArray object containing random values.

...
# Arguments
- `rown::Int64`: The number of rows in the array.
- `coln::Int64`: The number of columns in the array.
- `meancoloccs::Int64`: The mean number of non-zero elements in the array.
- `rowprefix::String`: The prefix to the row number. "Releve" by default.
- `colprefix::String`: The prefix to the column number. "Species" by default.
- `rowdim::String`: The row dimension name. "Releve" by default.
- `coldim::String`: The column dimension name. "Species" by default.
...

# Examples
```julia
julia>generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "Releve", colprefix = "Species")
```
"""
function generate_test_array(;rown::Int64, coln::Int64,
                             meancoloccs::Int64,
                             rowprefix::String = "Releve", colprefix::String = "Species",
                             rowdim::String = "Releve", coldim::String = "Species")

    # rown = 15
    # coln = 10
    # meancoloccs = 7
    # rowprefix = "SiteA-"
    # colprefix = "Species"
    # rowdim = "Releve"
    # coldim = "Species"

    zerop = meancoloccs / coln
    rownames = vec([string("$rowprefix")].*string.([1:1:rown;]))
    colnames = vec([string("$colprefix")].*string.([1:1:coln;]))
    x = NamedArrays.NamedArray(Array(sprand(Float64, rown, coln, zerop)), names = (rownames, colnames), dimnames = (rowdim, coldim))

    # x = NamedArrays.NamedArray(rand(2,2))
    # # x = rand(2, 2)
    # typeof(x)

    y = x ./ sum(x, dims = 2)

    return y

end

"""
Draws heavily from the function outlines here: https://discourse.julialang.org/t/nanmean-options/4994/17
"""
function nzfunc(f::Function, x::NamedArray; dims::Int64 = 1)

    # _nzfunc(fn, A, ::Colon) = fn(skip(iszero, A))
    # _nzfunc(fn, A) = map(a->_nzfunc(fn,a,:), eachcol(A))

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