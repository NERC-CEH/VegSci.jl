using NamedArrays
using SparseArrays
using InvertedIndices
using DataFrames
using Random

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
- `rowprefix::String`: The prefix to the row number. "Species" by default.
- `colprefix::String`: The prefix to the column number. "Releve" by default.
- `rowdim::String`: The row dimension name. "Species" by default.
- `coldim::String`: The column dimension name. "Releve" by default.
- `sparse_array::Bool`: If true a names sparse matrix is returned. If false a named dense matrix is returned. false by default.
...

# Examples
```julia
julia>VegSci.generate_test_array(rown = 10, coln = 10, meancoloccs = 5)
```
"""
function generate_test_array(;rown::Int64, coln::Int64, meancoloccs::Int64, 
                             val_type::String = "integer", scale::Union{Nothing, String} = "cols",
                             rowprefix::String = "Species", colprefix::String = "Releve",
                             rowdim::String = "Species", coldim::String = "Releve",
                             sparse_array::Bool = false)

    nonzerop = meancoloccs / coln
    rownames = vec([string("$rowprefix")].*string.([1:1:rown;]))
    colnames = vec([string("$colprefix")].*string.([1:1:coln;]))

    if val_type == "integer"
        digits = 0
        min_val = 1
        step_val = 1
        max_val = 100
        val_range = min_val:step_val:max_val
    elseif val_type == "float"
        digits = 2
        min_val = 0.01
        step_val = 0.01
        max_val = 1
        val_range = min_val:step_val:max_val
    elseif val_type == "bool"
        min_val = 0
        max_val = 1
        val_range = [min_val, max_val]
    end

    vals = sprand(rown, 1, 5 / rown, n -> rand(val_range, n))

    while nnz(vals) == 0
        vals = sprand(rown, 1, 5 / rown, n -> rand(val_range, n))
    end 

    for i in 2:coln

        vals_i = sprand(rown, 1, nonzerop, n -> rand(val_range, n))

        while nnz(vals_i) == 0
            vals_i = sprand(rown, 1, nonzerop, n -> rand(val_range, n))
        end 

        vals = hcat(vals, vals_i)

    end

    if sparse_array == true
        x = NamedArray(vals, names = (rownames, colnames), dimnames = (rowdim, coldim))
    elseif sparse_array == false
        x = NamedArray(Matrix(vals), names = (rownames, colnames), dimnames = (rowdim, coldim))
    end

    if val_type != "bool"

        if scale == "cols"
            scaling_vals = sum(x, dims = 1) * (1 / max_val)
        elseif scale == "rows"
            scaling_vals = sum(x, dims = 2) * (1 / max_val)
        end

        if sparse_array == true
            scaling_vals = sparse(scaling_vals)
        end

        if scale == "cols"
            scaling_mat = NamedArray(sparse(scaling_vals), names = (["sum"], colnames), dimnames = (["sum"], coldim))
        elseif scale == "rows"
            scaling_mat = NamedArray(sparse(scaling_vals), names = (rownames, ["sum"]), dimnames = (rowdim, ["sum"]))
        end

        x = x ./ scaling_mat

        x = round.(x, digits = digits)

        if val_type == "integer"
            x = Int.(x)
        end
        
    end

    return x

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