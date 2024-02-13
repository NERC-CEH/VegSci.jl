using NamedArrays
using SparseArrays
using InvertedIndices

"""
    greet_your_package_name(name)

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

    zerop = meancoloccs / coln
    rownames = vec([string("$rowprefix")].*string.([1:1:rown;]))
    colnames = vec([string("$colprefix")].*string.([1:1:coln;]))
    x = NamedArrays.NamedArray(Array(sprand(Float64, rown, coln, zerop)), names = (rownames, colnames), dimnames = (rowdim, coldim))
    y = x ./ sum(x, dims = 2)
    return y

end

"""
Draws heavily from the function outlines here: https://discourse.julialang.org/t/nanmean-options/4994/17
"""
function nzfunc(f::Function, x::NamedArray; dims::Int64 = 1)

    # _nzfunc(fn, A, ::Colon) = fn(skip(iszero, A))
    # _nzfunc(fn, A) = map(a->_nzfunc(fn,a,:), eachcol(A))
    # y = nzfunc(f, x, dims = dims)
    # setnames!(y, names(x)[2], 2)

    _nzfunc(fn, A, ::Colon) = fn(filter(!iszero, A))
    _nzfunc(fn, A, dims) = mapslices(a->_nzfunc(fn,a,:), A, dims=dims)
    nzfunc(fn, A; dims=:) = _nzfunc(fn, A, dims)
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