function indval_fidelity(x::NamedMatrix, clusters::Dict{String, Vector{String}})

    N = size(x)[1]
    n = sum(x, dims = 1)

    # if ismissing(clusters)

    #     Np = length(getindex(clusters, i))
    #     np = sum(x[getindex(clusters, i),:], dims = 1)
    #     indval = ((np .* (N - Np)) ./ (((n .* Np) .- (2 .* np)) .+ (np .* N))) .* (np ./ Np)

    # end

    indval_all = NamedArrays.NamedArray(zeros(length(names(clusters)), size(x)[2]), names = (names(clusters), names(x)[2]))

    for i in names(clusters)
        Np = length(getindex(clusters, i))
        np = sum(x[getindex(clusters, i),:], dims = 1)
        indval = ((np .* (N - Np)) ./ (((n .* Np) .- (2 .* np)) .+ (np .* N))) .* (np ./ Np)
        indval_all[i,:] = indval
    end

    return indval_all

end