using NamedArrays
using LinearAlgebra

"""
"""
function jaccard_similarity(x::Union{Matrix, NamedMatrix})
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
binary_similarity(x::NamedArray, eq::AbstractString)

Create a personalised greeting for VegSci using a `name`.

...
# Arguments
- `x::NamedArray`: A site by species 
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


function czekanowski_index(x::NamedArray, y::NamedArray)

    # Create a vector of pairwise samples and references to iterate over
    comp_vec = vec(collect(Iterators.product(names(x)[1], names(y)[1])))

    # Create an empty matrix to store results
    results = NamedArray(zeros(size(x, 1), size(y, 1)), names = (names(x)[1], names(y)[1]))

    # Loop through each pair of samples and references, calculate the Czekanowski index, and store the results in the matrix
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