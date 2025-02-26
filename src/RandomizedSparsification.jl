module RandomizedSparsification

using ProgressLogging, Logging
using LinearAlgebra, StatsBase, SparseArrays

include("utils.jl")
include("indices.jl")

@inline function _sparsify_single_threaded(
    A, A_1_norm, A_F_norm, P_ij, P_ij_cs; r, progress
) # *
    d₁, d₂ = size(A)
    probs = rand(Float32, r)
    ij = FindIndices._indices_groupsort(P_ij_cs, probs, d₁, d₂)
    vals = @inbounds (1 / r) * A[ij] ./ P_ij[ij]
    A_hat = sparse(ij, vals, d₁, d₂)
    return A_hat
end

@inline function _sparsify_multi_threaded(A, A_1_norm, A_F_norm, P_ij, P_ij_cs; r, progress) # *
    d₁, d₂ = size(A)
    counts = Dict{CartesianIndex,Int32}()
    for _ in Base.OneTo(r)
        c = rand()
        ind = CartesianIndices(A)[searchsortedfirst(P_ij_cs, c)] # **
        if !haskey(counts, ind)
            counts[ind] = 1
        else
            counts[ind] += 1
        end
    end
    # counts = StatsBase.countmap(ij)
    vals = map([keys(counts)...]) do ij
        @inbounds (counts[ij] / r) * A[ij] / P_ij[ij] # ***
    end
    A_hat = sparse([keys(counts)...], vals, d₁, d₂) # ****
    return A_hat
end

"""
	_sparsify_context_data(A)

Compute the necessary context data, `A_1_norm`, `A_F_norm`, `P_ij` and `P_ij_cs` (cumulative distribution on the matrix elements) for the `sparsify` randomized sparisification algorithm of `A`.
"""
@inline function _sparsify_context_data(A)
    A_1_norm = norm(A, 1)
    A_F_norm = norm(A, 2) # Frobenius norm is equivalent to the p=2 norm
    P_ij = @. 1 / 2 * (abs(A)^2 / A_F_norm^2 + abs(A) / A_1_norm) # Probabilities of the elements
    P_ij_cs = cumsum(view(P_ij, :)) # *****
    return A_1_norm, A_F_norm, P_ij, P_ij_cs
end

"""
	sparsify(A; r, multithreaded=false, progress=false)

Sparsify the input matrix `A` using the randomized element-wise sparsification algorithm with `r` samples.
"""
function sparsify(A; r, multithreaded=false, progress=false)
    ctxt_data = _sparsify_context_data(A)
    return if multithreaded
        _sparsify_multi_threaded(A, ctxt_data...; r, progress)
    else
        _sparsify_single_threaded(A, ctxt_data...; r, progress)
    end
end

export sparsify
end
