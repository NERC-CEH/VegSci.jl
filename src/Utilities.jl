using NamedArrays

"""
    greet_your_package_name(name)

Create a `rown` by `coln` NamedArray object containing random values.

...
# Arguments
- `rown::Int64`: The number of rows in the array.
- `coln::Int64`: The number of columns in the array.
- `min::Float64`: The minimum possible random value. 0.0 by default.
- `max::Float64`: The maximum possible random value. 1.0 by default.
- `increment::Float64`: The increment used in generating random values. 0.1 by default.
- `rowprefix::String`: The prefix to the row number. "Releve" by default.
- `colprefix::String`: The prefix to the column number. "Species" by default.
...

# Examples
```julia
julia>generate_test_array(rown = 10, coln = 10, min = 0.0, max = 1.0, increment = 0.1, rowprefix = "Releve", colprefix = "Species")
```
"""
function generate_test_array(;rown::Int64, coln::Int64, 
                             min::Float64 = 0.0, max::Float64 = 1.0, increment::Float64 = 0.1, 
                             rowprefix::String = "Releve", colprefix::String = "Species")

    rownames = vec([string("$rowprefix")].*string.([1:1:rown;]))
    colnames = vec([string("$colprefix")].*string.([1:1:coln;]))
    x = NamedArrays.NamedArray(rand(min:increment:max, rown, coln), names = (rownames, colnames), dimnames = (rowprefix, colprefix))
    return x

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