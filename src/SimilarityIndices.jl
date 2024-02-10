using NamedArrays
using LinearAlgebra

"""
    greet_your_package_name(name)

Create a personalised greeting for EcoVeg using a `name`.

...
# Arguments
- `x::NamedArray`: A site by species 
...

# Examples
```julia
x = generatetestarray(rown = 10, coln = 10)
julia>binarysimilarity(x, "(a ./ (a .+ b .+ c)) + I")
```
"""
function binary_similarity(x::NamedArray, eq::AbstractString = "(a ./ (a .+ b .+ c)) + I")

	samp = Int.(x .!= 0)
	d = samp * transpose(samp)
	s = Array(diag(d))
	n = length(s)
	a = d - Diagonal(d)
	b = reshape(repeat(s, n), :, n) - a
	c = Array(transpose(reshape(repeat(s, n), :, n))) - a

	# sim = eval(Meta.parse(eq))

    expr = Meta.parse(eq)

    sim = eval(:((a,b,c)->$expr))

	return sim

end

function czekanowski_index(x::NamedArray, y::NamedArray)

    # Create a vector of pairwise samples and references to iterate over
    comp_vec = vec(collect(Iterators.product(names(x)[1], names(y)[1])))

    # Create an empty matrix to store results
    results = NamedArray(zeros(size(x, 1), size(y, 1)), names = (names(x)[1], names(y)[1]))

    # Loop through each pair of samples and references, calculate the Czekanowski index, and store the results in the matrix
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