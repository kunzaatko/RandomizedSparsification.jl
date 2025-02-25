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
