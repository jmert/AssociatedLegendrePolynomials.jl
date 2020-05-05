"""
    abstract type AbstractLegendreNorm end

Abstract trait supertype for normalization conditions of the Associated Legendre polynomials.
"""
abstract type AbstractLegendreNorm end

# Define element types as a way to conveniently promote, even for the trait types.
# Union{} will always lose out in type promotion.
Base.eltype(::Type{<:AbstractLegendreNorm}) = Union{}
# Normalizations are scalars in a broadcasting context
Base.broadcastable(x::AbstractLegendreNorm) = Ref(x)

## Below are overloadable non-exported functions to implement the AbstractLegendreNorm
## interface.

"""
    Plm_00(::N, ::Type{T}) where {N<:AbstractLegendreNorm, T}

Returns the initial condition ``P_0^0(x)`` for the associated Legendre recursions based
on the normalization choice `N` for numeric type `T`.
"""
function Plm_00 end

"""
    Plm_μ(norm::N, ::Type{T}, l::Integer) where {N<:AbstractLegendreNorm, T}

Returns the coefficient ``μ_ℓ`` for the single-term recursion relation
```math
    P_{ℓ+1}^{ℓ+1}(x) = -μ_{ℓ+1} \\sqrt{1-x^2} P_ℓ^ℓ(x)
```
where ``μ_ℓ`` is appropriate for the choice of normalization `N`.
"""
function Plm_μ end

"""
    Plm_ν(norm::N, ::Type{T}, l::Integer) where {N<:AbstractLegendreNorm, T}

Returns the coefficient ``ν_ℓ`` for the single-term recursion relation
```math
    P_{ℓ+1}^ℓ(x) = ν_{ℓ+1} x P_ℓ^ℓ(x)
```
where ``ν_ℓ`` is appropriate for the choice of normalization `N`.
"""
function Plm_ν end

"""
    Plm_α(norm::N, ::Type{T}, l::Integer, m::Integer) where {N<:AbstractLegendreNorm, T}

Returns the coefficient ``α_ℓ^m`` for the two-term recursion relation
```math
    P_{ℓ+1}^{m}(x) = α_{ℓ+1}^m x P_ℓ^m(x) - β_{ℓ+1}^m P_{ℓ-1}^m(x)
```
where ``α_ℓ^m`` is appropriate for the choice of normalization `N`.
"""
function Plm_α end

"""
    Plm_β(norm::N, ::Type{T}, l::Integer, m::Integer) where {N<:AbstractLegendreNorm, T}

Returns the coefficient ``β_ℓ^m`` for the two-term recursion relation
```math
    P_{ℓ+1}^{m}(x) = α_{ℓ+1}^m x P_ℓ^m(x) - β_{ℓ+1}^m P_{ℓ-1}^m(x)
```
where ``β_ℓ^m`` is appropriate for the choice of normalization `N`.
"""
function Plm_β end

"""
    boundscheck_hook(norm::AbstractLegendreNorm, lmax, mmax)

A bounds-checking hook executed at the beginning of each [`legendre!`](@ref) call to
permit a normalization `norm` to validate that the given maximum ``(ℓ,m)`` will be within
the ability to satisfy. The default case always returns `nothing`. A custom normalization
should throw an error if `lmax` or `mmax` is out of bounds or return `nothing` otherwise.

For example, the precomputed coefficients of [`LegendreNormCoeff`](@ref) are limited to
a given domain at time of construction and cannot be used to calculate terms to arbitrary
orders/degrees.
"""
function boundscheck_hook end

boundscheck_hook(norm::AbstractLegendreNorm, lmax, mmax) = nothing
