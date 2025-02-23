module RandomizedSparsification

using ProgressLogging, Logging
using LinearAlgebra, StatsBase, SparseArrays

import SparseArrays: sparse

# NOTE: Adding upstream: https://github.com/JuliaSparse/SparseArrays.jl/pull/600 <22-02-25> 
# FIX: Known type piracy. When added upstream, this method can be removed. <22-02-25> 
function sparse(IJ::Vector{<:CartesianIndex}, v, m, n)
    IJ′ = reinterpret(Int, reshape(IJ, 1, :))
    return sparse(view(IJ′, 1, :), view(IJ′, 2, :), v, m, n)
end

"""
	rankrand(d₁, d₂, rank)
Generate a random matrix of dimensions `d₁` × `d₂` with a rank `rank`.
"""
function rankrand(d₁, d₂, rank)
    @assert max(d₁, d₂) >= rank "The rank will not be correct because the specified dimension is not sufficient"
    basevecs = randn(max(d₁, d₂), rank)
    A = basevecs * basevecs'
    return A[1:d₁, 1:d₂]
end

@inline function _sparsify_single_threaded(
    A, A_1_norm, A_F_norm, P_ij, P_ij_cs; r, progress
) # *
    d₁, d₂ = size(A)
    counts = Dict{CartesianIndex,Int32}()
    probs = rand(Float32, r)
    Logging.with_logger(progress ? Logging.current_logger() : Logging.NullLogger()) do
        @progress for c in probs
            ind = CartesianIndices(A)[searchsortedfirst(P_ij_cs, c)] # **
            if !haskey(counts, ind)
                counts[ind] = 1
            else
                counts[ind] += 1
            end
        end
    end
    # counts = StatsBase.countmap(ij)
    vals = map([keys(counts)...]) do ij
        @inbounds (counts[ij] / r) * A[ij] / P_ij[ij] # ***
    end
    A_hat = sparse([keys(counts)...], vals, d₁, d₂) # ****
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

@doc """
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

@inline r_ε1(d₁, d₂, ε, rank) = ε^(-2) * max(d₁, d₂) * log(d₁ + d₂) * rank
@inline r_ε2(d₁, d₂, ε, rank) = 2 / 3 * ε^(-1) * log(d₁ + d₂) * √(d₁ * d₂) * √(rank)

@doc raw"""
	srank(A::AbstractMatrix)

Compute the stable rank of the input matrix `A`.

The stable rank is defined as the ratio of the squared Frobenius norm to the squared spectral norm of the matrix.
```math
	\frac{\|A\|^2_\text{F}}{\|A\|^2_\text{op}}
```
"""
srank(A::AbstractMatrix) = norm(A, 2)^2 / opnorm(A, 2)^2

function r_ε(A::AbstractMatrix, ε::Real, rank=nothing)
    return r_ε(size(A)..., ε, !isnothing(rank) ? rank : srank(A))
end

"""
	r_ε(A, ε, [rank])
	r_ε(d₁, d₂, ε, rank)

Compute the sufficient number of samples `r` required to achieve a given error tolerance `ε` in the sparsification of the matrix `A` (``A \\in \\mathbb{R}^{d_1 \\times d_2}``) in the expected value. Optionally use `rank` as an upperbound for ``\\operatorname{srank}``. The number of samples returned by this function is not necessarily the least bound for the number of samples needed for the error to bound to hold.
"""
function r_ε(d₁::Integer, d₂::Integer, ε::Real, rank)
    return Int(ceil(max(r_ε1(d1, d₂, ε, rank), r_ε2(d1, d₂, ε, rank))))
end

@doc raw"""
	ρ_s(A::AbstractSparseArray)

Compute the proportional sparsity of the input matrix `A`.
```math
	\frac{\#\{(i,j)\mid A_{ij} \neq 0 \}}{\#\{(i,j)\}}
```
"""
ρ_s(A::AbstractSparseArray) = length(A.nzval) / prod(size(A))

@doc raw"""
	representation_error(A, A_r)

Compute the relative operator representation error between the original matrix `A` and the sparsified matrix `A_r`.
```math
	\frac{\|A - A_r\|_\text{op}}{\|A\|_\text{op}}
```
"""
representation_error(A, A_r) = opnorm(A - A_r) / opnorm(A)
end
