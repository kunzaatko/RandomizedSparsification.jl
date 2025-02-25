"""
This module contains functions for determining the indices from the cumulative sum of the probability matrix and the
randomly sampled numbers from ``[0,1]`` of the length of the desired number of samples.

It does not consider whether or not the number of the samples that are taken fit into the memory or not. This should be
considered elsewhere in the call stack where the number of samples desired can be separated into batches that should be
sampled at once and only the matrix should be held between the iterations.

The functions that find the summation indices in the sparsification algorithm have the signature signature: 
```julia
_indices(
    P_ij_cs::AbstractVector{<:Real}, 
    rand01_vec::AbstractVector{<:Real},
    m::Integer, n::Integer
)::AbstractDict{<:CartesianIndex, <:Unsigned}
```
where `P_ij_cs` is the vector of the cumulative sum of the probabilities, `rand01_vec` are the sampled values and `m`
and `n` are the dimensions of the original matrix. So the indices are returned in the form of
a `<:AbstractDict{<:CartesianIndex, <:Unsigned}`, where the `Unsigned` holds the counts of the samples of
`ind::CartesianIndex`.
""" 
module FindIndices
@inline function _indices_sequencial()
end

@inline function _indices_pivoting()
end
end
