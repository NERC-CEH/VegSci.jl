function check_cluster_releveID_duplicates(clusters::Dict)

    clusters = copy(clusters)

    duplicates = DataFrames.DataFrame(releveSamp = [], releveComp = [], duplicates = [])

    for i in names(clusters)

        i_names = getindex(clusters, i)

        for j in setdiff(names(clusters), [i])

        j_names = getindex(clusters, j)

        releve_duplicates = collect(intersect(Set(i_names), Set(j_names)))

        if length(releve_duplicates) > 0

            for k in releve_duplicates

                duplicates_ijk = DataFrames.DataFrame(releveSamp = i, releveComp = j, duplicates = k)

                append!(duplicates, duplicates_ijk)

            end

        end

        delete!(clusters, i)

        end

    end

    return duplicates

end

"""
    check_relSpp_array_format(x::NamedArray; lower_bound::Float64 = 0.0, upper_bound::Float64 = 1.0)

Check whether the format of a Releve by Species matrix is correct. 
Specifically whether the following conditions are met:
    -   All column names are strings
    -   All row names are strings
    -   Cover-Abundance values are above or equal to a lower bound
    -   Cover-Abundance values are below or equal to a upper bound
    -   There are no NaN Cover-Abundance values
    -   There are no missing Cover-Abundance values

...
# Arguments
- `x::NamedArray`: A Releve by Species matrix
- `lower_bound::Float64`: The minimum allowable cover-abundance value.
- `upper_bound::Float64`: The maximum allowable cover-abundance value.
...

# Examples
```julia
julia> x = generate_test_array(rown = 10, coln = 10, meancoloccs = 5, rowprefix = "Releve", colprefix = "Species")
julia> x = check_relSpp_array_format(x)
```
"""
function check_relSpp_array_format(x::NamedArray; lower_bound::Float64 = 0.0, upper_bound::Float64 = 1.0)

    # Check that column names are all of type String
    column_names_strings = typeof(names(x, 2)) <: Vector{String}

    # Check that row names are all of type String
    row_names_strings = typeof(names(x, 1)) <: Vector{String}

    # Check that cover-abundance values are between a lower and upper bounds
    coverabundance_above_lowerbound = all(x .>= lower_bound)
    coverabundance_below_upperbound  = all(x .<= upper_bound)

    # Check that cover-abundance values do not contain NaN
    coverabundance_no_NaN = !any(isnan, x)

    # Check that cover-abundance values do not contain missing
    coverabundance_no_missing = !any(ismissing, x)

    # Check whether all conditions are true
    correct_format = all([column_names_strings,
                          row_names_strings,
                          coverabundance_above_lowerbound,
                          coverabundance_below_upperbound,
                          coverabundance_no_NaN,
                          coverabundance_no_missing
                          ])

    # Compile checks
    checks = Dict("column_names_strings" => column_names_strings,
                  "row_names_strings" => row_names_strings,
                  "coverabundance_above_lowerbound" => coverabundance_above_lowerbound,
                  "coverabundance_below_upperbound" => coverabundance_below_upperbound,
                  "coverabundance_no_NaN" => coverabundance_no_NaN,
                  "coverabundance_no_missing" => coverabundance_no_missing,
                  "correct_format" => correct_format
                  )

    return checks

end