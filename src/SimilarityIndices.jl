using NamedArrays
using LinearAlgebra
using SparseArrays
using VegSci
"""
"""
function jaccard_coefficient(x::Union{Matrix, NamedMatrix})
    samp = Int.(x .!= 0)
	d = samp * transpose(samp)
	s = Array(diag(d))
	n = length(s)
	a = d - Diagonal(d)
	b = reshape(repeat(s, n), :, n) - a
	c = Array(transpose(reshape(repeat(s, n), :, n))) - a
	j = (a ./ (a .+ b .+ c)) + I
	return j    
end

"""
binary_dissimilarity(x::Union{Matrix, NamedMatrix}, eq::AbstractString)

Create a personalised greeting for VegSci using a `name`.

...
# Arguments
- `x::Union{Matrix, NamedMatrix}`: A releve by species matrix of the class Matrix or NamedArrays.NamedMatrix
- `eq::AbstractString`: A representing an equation using the terms a, b, and c. The Jaccard similarity by default. "(a ./ (a .+ b .+ c)) + I".
...

# Examples
```julia
x = generatetestarray(rown = 10, coln = 10)
julia>binarysimilarity(x, "(a ./ (a .+ b .+ c)) + I")
```
"""
function binary_dissimilarity(x::Union{Matrix, NamedMatrix}, eq::String)
	samp = Int.(x .!= 0)
	d = samp * transpose(samp)
	s = Array(diag(d))
	n = length(s)
	a = d - Diagonal(d)
	b = reshape(repeat(s, n), :, n) - a
	c = Array(transpose(reshape(repeat(s, n), :, n))) - a
	sim = eval(Meta.parse(eq))

	f(a,b,c) = eval(Meta.parse(eq))
	sim = f(a,b,c)
	return sim
end


function steinhaus_coefficient(x::NamedArray, y::NamedArray)

    # Create a vector of pairwise samples and references to iterate over
    comp_vec = vec(collect(Iterators.product(names(x)[1], names(y)[1])))

    # Create an empty matrix to store results
    results = NamedArray(zeros(size(x, 1), size(y, 1)), names = (names(x)[1], names(y)[1]))

    # Loop through each pair of samples and references, calculate the Steinhaus coefficient, and store the results in the matrix
	# There will be a much better way to iterate through each pair of samples and references!
    for i in comp_vec
        x_i = x[[i[1]],:]
        y_i = y[[i[2]],:]
        eval_mat = vcat(x_i, y_i)
        A = sum(x_i, dims = 2)[1]
        B = sum(y_i, dims = 2)[1]
        W = sum(minimum(eval_mat, dims = 1), dims = 2)[1,1]
        sim = 2 * W / (A + B)
        results[i[1], i[2]] = sim 
    end

    return results

end


# """
# similarity_fqi(x::NamedMatrix)::NamedMatrix

# Calculate the Frequency Index (FQI) between two...

# ### Input

# - `x` -- A releve by species matrix, of class NamedArrays::NamedMatrix
# - `y` -- A vegetation unit by species matrix, of class NamedArrays::NamedMatrix

# ### Output

# The FQI, a value between 0 and 100.

# ### Notes


# ### Algorithm

# This function implements 

# ### References

# Tichý, L., 2005. New similarity indices for the assignment of relevés to the vegetation units of an existing phytosociological classification. Plant Ecol 179, 67–72. https://doi.org/10.1007/s11258-004-5798-8

# """
# function similarity_fqi(x::NamedMatrix, y::NamedMatrix)
# 	x = VegSci.generate_test_array(rown = 10, coln = 20, meancoloccs = 15, rowprefix = "SiteX-", colprefix = "Species")[:,Not(["Species5", "Species6", "Species20"])]
# 	y = VegSci.generate_test_array(rown = 5, coln = 20, meancoloccs = 15, rowprefix = "SiteY-", colprefix = "Species")[:,Not(["Species3", "Species11", "Species17"])]
# 	cols = collect(intersect(Set(names(y)[2]), Set(names(x)[2])))
# 	x_sum = sum(y[:, cols], dims = 2)
# 	y_sum = sum(y, dims = 2)
# 	fqi = 100 .* (x_sum ./ y_sum)
# 	setnames!(fqi, names(x)[1], 2)
# 	setdimnames!(fqi, ["y", "x"])

# 	return fqi

# end

# function similarity_pfdi(x::NamedMatrix, y::NamedMatrix)

# 	cols = collect(intersect(Set(names(y)[2]), Set(names(x)[2])))
# 	x_sum = sum(y[:, cols], dims = 2)
# 	y_sum = sum(y, dims = 2)
# 	fqi = 100 .* (x_sum ./ y_sum)
# 	setnames!(fqi, names(x)[1], 2)
# 	setdimnames!(fqi, ["y", "x"])

# 	return fqi

# end