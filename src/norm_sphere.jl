# Implements the legendre interface for the unit normalization case

"""
    struct LegendreOrthoNorm <: AbstractLegendreNorm end

Trait type denoting the orthonormal (full) normalization of the associated Legendre
functions ``O_\\ell^m(x)``.
The orthonormal normalization is defined such that
```math
    \\int_{-1}^{1} \\left[ O_\\ell^m(x) \\right]^2 \\,dx = 1
```
The normalization factor with respect to the standard (unnormalized) Legendre
functions ``P_\\ell^m(x)`` ([`LegendreUnitNorm`](@ref)) is given by
```math
    O_\\ell^m(x) \\equiv \\sqrt{\\frac{2\\ell+1}{2} \\frac{(\\ell-m)!}{(\\ell+m)!}} P_\\ell^m(x)
```
"""
struct LegendreOrthoNorm <: AbstractLegendreNorm end

"""
    struct LegendreSphereNorm <: AbstractLegendreNorm end

Trait type denoting the spherical-harmonic normalization of the associated Legendre
functions ``\\lambda_\\ell^m(x)``.
The spherical-harmonic normalization is defined such that
```math
    \\int_{-1}^{1} \\left[ \\lambda_\\ell^m(x) \\right]^2 \\,dx = \\frac{1}{2\\pi}
```
The normalization factor with respect to the standard (unnormalized) Legendre
functions ``P_\\ell^m(x)`` ([`LegendreUnitNorm`](@ref)) is given by
```math
    \\lambda_\\ell^m(x) \\equiv \\sqrt{\\frac{2\\ell+1}{4\\pi} \\frac{(\\ell-m)!}{(\\ell+m)!}} P_\\ell^m(x)
```
"""
struct LegendreSphereNorm <: AbstractLegendreNorm end

@inline function initcond(::LegendreOrthoNorm, ::Type{T}) where T
    return sqrt(inv(T(2)))
end
@inline function initcond(::LegendreSphereNorm, ::Type{T}) where T
    # comparing this against
    #   convert(T, inv(sqrt(4*convert(BigFloat, π))))
    # shows that this is exact within Float64 precision
    return inv(sqrt(4 * convert(T, π)))
end

const OrthoOrSphereNorm = Union{LegendreOrthoNorm,LegendreSphereNorm}

# Version of sqrt() which skips the domain (x < 0) check for the IEEE floating point types.
# For nonstandard number types, just fall back to a regular sqrt() since eliminating the
# domain check is probably no longer the dominant contributor to not vectorizing.
unchecked_sqrt(x::T) where {T <: Base.IEEEFloat} = Base.sqrt_llvm(x)
unchecked_sqrt(x::T) where {T <: Integer} = unchecked_sqrt(float(x))
unchecked_sqrt(x) = Base.sqrt(x)

@inline function
coeff_μ(::OrthoOrSphereNorm, ::Type{T}, l::Integer) where T
    # The direct derivation of the normalization constant gives
    #     return sqrt(one(T) + inv(convert(T, 2l)))
    # but when comparing results for T ∈ (Float64,BigFloat), the Float64 results differ by
    # ±1 ulp at various points.
    #
    # Instead, consider the following which limits to 0 (i.e. as (2l)^-1) rather than 1.
    # ```math
    #     μ'_ℓ = μ_ℓ - 1 = \sqrt{1 + \frac{1}{2ℓ}} - 1
    # ```
    # Now complete the square by multiplying through by a factor of unity to instead write
    # ```math
    #   \left( \sqrt{1 + \frac{1}{2ℓ}} - 1 \right) \times
    #       \frac{\sqrt{1 + \frac{1}{2ℓ}} + 1}{\sqrt{1 + \frac{1}{2ℓ}} + 1}
    #   = \frac{1}{\sqrt{2ℓ(2ℓ + 1) + 2ℓ}}
    # ```
    # Then rewrite ``μ_ℓ`` in terms of ``μ'_ℓ`` to arrive at
    # ```math
    #   μ_ℓ = 1 + \left[ 2ℓ + \sqrt{2ℓ(2ℓ + 1)} \right]^{-1}
    # ```
    # The square root is calculated on a growing quantity rather than one limiting to 1,
    # so the significance of the low-order digits are better preserved. This expression
    # allows Float64 calculations to be correctly rounded values when compared to BigFloat
    # to very high ℓ.
    xT = convert(T, 2l)
    return one(T) + inv(xT + unchecked_sqrt(muladd(xT, xT, xT)))
end

@inline function
coeff_ν(::OrthoOrSphereNorm, ::Type{T}, l::Integer) where T
    return unchecked_sqrt(convert(T, 2l + 1))
end

@inline function
coeff_α(::OrthoOrSphereNorm, ::Type{T}, l::Integer, m::Integer) where T
    lT = convert(T, l)
    mT = convert(T, m)
    # Write factors in two pieces to make compiler's job easier. In the case where
    # both coeff_α and coeff_β are called and inlined, the next line from both functions
    # should be merged and shared.
    fac1 = (2lT + 1) / ((2lT - 3) * (lT^2 - mT^2))
    fac2 = 4*(lT - 1)^2 - 1
    return unchecked_sqrt(fac1 * fac2)
end

@inline function
coeff_β(::OrthoOrSphereNorm, ::Type{T}, l::Integer, m::Integer) where T
    lT = convert(T, l)
    mT = convert(T, m)
    # Write factors in two pieces to make compiler's job easier. In the case where
    # both coeff_α and coeff_β are called and inlined, the next line from both functions
    # should be merged and shared.
    fac1 = (2lT + 1) / ((2lT - 3) * (lT^2 - mT^2))
    fac2 = (lT - 1)^2 - mT^2
    return unchecked_sqrt(fac1 * fac2)
end

# Extra functions

"""
    N = Nlm([T=Float64], l, m)

Computes the normalization constant
```math
    N_\\ell^m \\equiv \\sqrt{\\frac{2\\ell+1}{4\\pi} \\frac{(\\ell-m)!}{(\\ell+m)!}}
```
which defines the Spherical Harmonic normalized functions ``\\lambda_\\ell^m(x)`` in
terms of the standard unit normalized ``P_\\ell^m(x)``
```math
    \\lambda_\\ell^m(x) \\equiv N_\\ell^m P_\\ell^m(x)
```
using numbers of type `T`.

See also [`Plm`](@ref) and [`λlm`](@ref).
"""
function Nlm(::Type{T}, l::Integer, m::Integer) where T
    fac1 = one(T)
    for ii in (l-m+1):(l+m)
        fac1 *= convert(T, ii)
    end
    num = 2*convert(T, l) + 1
    den = 4*convert(T, π)
    return sqrt( num * inv(den) * inv(fac1) )
end
Nlm(l::Integer, m::Integer) = Nlm(Float64, l, m)
