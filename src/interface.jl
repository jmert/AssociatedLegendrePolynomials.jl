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
    initcond(::N, ::Type{T}) where {N<:AbstractLegendreNorm, T}

Returns the initial condition ``P_0^0(x)`` for the associated Legendre recursions based
on the normalization choice `N` for numeric type `T`.
"""
function initcond end

"""
    coeff_μ(norm::N, ::Type{T}, l::Integer) where {N<:AbstractLegendreNorm, T}

Returns the coefficient ``\\mu_\\ell`` for the single-term recursion relation
```math
    P_\\ell^\\ell(x) = -\\mu_\\ell \\sqrt{1-x^2} P_{\\ell-1}^{\\ell-1}(x)
```
where ``\\mu_\\ell`` is appropriate for the choice of normalization `N`.
"""
function coeff_μ end

"""
    coeff_ν(norm::N, ::Type{T}, l::Integer) where {N<:AbstractLegendreNorm, T}

Returns the coefficient ``\\nu_\\ell`` for the single-term recursion relation
```math
    P_\\ell^{\\ell-1}(x) = \\nu_\\ell x P_{\\ell-1}^{\\ell-1}(x)
```
where ``\\nu_\\ell`` is appropriate for the choice of normalization `N`.
"""
function coeff_ν end

"""
    coeff_α(norm::N, ::Type{T}, l::Integer, m::Integer) where {N<:AbstractLegendreNorm, T}

Returns the coefficient ``\\alpha_\\ell^m`` for the two-term recursion relation
```math
    P_\\ell^m(x) = \\alpha_\\ell^m x P_{\\ell-1}^m(x) - \\beta_\\ell^m P_{\\ell-2}^m(x)
```
where ``\\alpha_\\ell^m`` is appropriate for the choice of normalization `N`.
"""
function coeff_α end

"""
    coeff_β(norm::N, ::Type{T}, l::Integer, m::Integer) where {N<:AbstractLegendreNorm, T}

Returns the coefficient ``\\beta_\\ell^m`` for the two-term recursion relation
```math
    P_\\ell^m(x) = \\alpha_\\ell^m x P_{\\ell-1}^m(x) - \\beta_\\ell^m P_{\\ell-2}^m(x)
```
where ``\\beta_\\ell^m`` is appropriate for the choice of normalization `N`.
"""
function coeff_β end

"""
    boundscheck_hook(norm::AbstractLegendreNorm, lmax, mmax)

A bounds-checking hook executed at the beginning of each [`legendre!`](@ref) call to
permit a normalization `norm` to validate that the given maximum ``(\\ell,m)`` will be
within the ability to satisfy. The default case always returns `nothing`. A custom
normalization should throw an error if `lmax` or `mmax` is out of bounds or return `nothing`
otherwise.

For example, the precomputed coefficients of [`LegendreNormCoeff`](@ref) are limited to
a given domain at time of construction and cannot be used to calculate terms to arbitrary
orders/degrees.
"""
function boundscheck_hook end

boundscheck_hook(norm::AbstractLegendreNorm, lmax, mmax) = nothing
