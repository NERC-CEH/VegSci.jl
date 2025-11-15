using NamedArrays

"""
    hamming_distance(nm::NamedMatrix)::NamedMatrix

Calculate the Hamming distance between all releves present in a species by releve matrix.

### Input

- `nm` -- A species by releve matrix containing categorical values, of class NamedArrays::NamedMatrix

### Output

A releve by releve NamedMatrix containing the Hamming Distance between each pair of releves.

### Notes

### Algorithm

This function calculates the Hamming distance as describes in Hamming (1950), adapted from https://johanndejong.wordpress.com/2019/09/28/fast-weighted-hamming-distance-in-r/.

### References

Hamming, R.W., 1950. Error Detecting and Error Correcting Codes. Bell System Technical Journal 29, 147-160. https://doi.org/10.1002/j.1538-7305.1950.tb00463.x

"""
function hamming_distance(nm::NamedMatrix)

    uniq_vals = unique(vec(nm))
    U = nm .== uniq_vals[1]
    H = U' * U

    for val in uniq_vals[2:end]
        U = nm .== val
        H = H + U' * U
    end

    hamming_dist = size(nm)[1] .- H

    return hamming_dist

end