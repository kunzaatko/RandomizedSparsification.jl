"""
This module contains functions for determining the indices from the cumulative sum of the probability matrix `P_ij_cs`
(i.e. cumulative distribution function on the indices) and the randomly sampled numbers from ``[0,1]`` of the length of
the desired number perspired by the user to align with the necessary accuracy bounds.

!!! danger
    The methods do not check whether or not the number of the samples that are taken fit into the RAM memory. This
    should be considered before in the call stack. The number of samples desired must be separated into batches that fit
    into the memory of the particular algorithm used for the indices location.

The functions that find the summation indices in the sparsification algorithm have the signature: 
```julia
_indices(
    P_ij_cs::AbstractVector{<:Real},    # cumulative distribution function on matrix indices
    rand01_vec::AbstractVector{<:Real}, # sampled values in [0,1]
    m::Integer, n::Integer              # original matrix size
)::AbstractVector{CartesianIndex{2}}
```
So the indices are returned as the vector `[ij::CartesianIndex{2}...]` with repetitions. This is passed to the
[`sparse`](@extref SparseArrays.sparse) function for the creation of the sparsified matrix.
"""
module FindIndices
@inline function _indices_sequencial() end

@inline function _indices_pivoting() end
end
