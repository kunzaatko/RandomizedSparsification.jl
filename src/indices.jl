"""
This module contains functions for determining the indices from the cumulative sum of the probability matrix `P_ij_cs`
(i.e. cumulative distribution function on the indices) and the randomly sampled numbers from ``[0,1]`` of the length of
the desired number perspired by the user to align with the necessary accuracy bounds.


!!! danger
    1. The methods do not check whether or not the number of the samples that are taken fit into the RAM memory.
        This should be considered before in the call stack. The number of samples desired must be separated into
        batches that fit into the memory of the particular algorithm used for the indices location.
    2. Assumption that `P_ij_cs` fits into memory is made. If it does not we need to be able to compute it at evaluation
        time for a given index efficiently (which is [not yet
        implemented](https://github.com/kunzaatko/RandomizedSparsification.jl/issues/1)).

The functions that find the summation indices in the sparsification algorithm have the signature: 
```julia
_indices(
    P_ij_cs::AbstractVector{<:Real},    # cumulative distribution function on matrix indices
    rands01::AbstractVector{<:Real},    # sampled values in [0,1]
    m::Integer, n::Integer;             # original matrix size
    progress=false                      # report progress to the user
)::AbstractVector{CartesianIndex{2}}
```
So the indices are returned as the vector `[ij::CartesianIndex{2}...]` with repetitions. This is passed to the
[`sparse`](@extref SparseArrays.sparse) function for the creation of the sparsified matrix.
"""
module FindIndices
using DataFrames
@inline function _indices_sequencial(P_ij_cs, rands01, m, n; progress=false)
    # FIX: suppresses all the logging. Instead I would like to suppress only the progress logging. Is there a way of
    # doing that using `ProgressLogging.jl`? <25-02-25> 
    # Logging.with_logger(progress ? Logging.current_logger() : Logging.NullLogger()) do
    # @progress for c in 
    return map(p -> CartesianIndices((m, n))[searchsortedfirst(P_ij_cs, p)], rands01)
end

"""
    _indices_groupsort(P_ij_cs, rands01, m, n)
    
# Algorithm
1. Defines a [`DataFrame`](@extref DataFrames.DataFrame) holiding `P_ij_cs` and `rands01` in the `:vals` column having a boolean `:group` column that differentiates values from `P_ij_cs` from `rands01`.
2. Rows are sorted by the `:vals` columns interleaving the `rands01`.
3. Indices are determined by mapping the number of preceding "`:group == P_ij_cs`" values for ta given `:group == rands01`

!!! perf 
        `DataFrames.jl` should only bring forward a small overhead compared to doing the underlying sorting manually.
"""
@inline function _indices_groupsort(P_ij_cs, rands01, m, n; progress=false) 
    df = DataFrame(:vals=>[P_ij_cs; rands01], :group=>[fill(true, length(P_ij_cs)); fill(false, length(rands01))])
    # TODO: Choose algorithm for sorting that uses the fact that the vector of P_ij_cs is already sorted. <26-02-25> 
    sort!(df, [:vals, :group])
    df[!, :ind] = cumsum(df[!,  :group])
    return getindex(CartesianIndices((m,n)), df[df[!, :group] .== false ,:ind] .+ 1)
end

# TODO: Recursive function that gives the index when it finds it otherwise it splits the input into two and passes the
# splits to the function once it otherwise it splits the input into two and passes the splits to the function once again <26-02-25> 
@inline function _indices_pivoting(P_ij_cs, rands01, m, n; progress=false) end

# TODO: Two indexed piping indexes algorithm <26-02-25> 
end
