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