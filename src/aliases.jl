"""
    LegendreUnitCoeff{T}

Precomputed recursion relation coefficients for the standard unit
normalization. Alias for `LegendreNormCoeff{LegendreUnitNorm,T}`.
"""
LegendreUnitCoeff{T} = LegendreNormCoeff{LegendreUnitNorm,T}

"""
    LegendreSphereCoeff{T}

Table type of precomputed recursion relation coefficients for the spherical
harmonic normalization. Alias for `LegendreNormCoeff{LegendreSphereNorm,T}`.
"""
LegendreSphereCoeff{T} = LegendreNormCoeff{LegendreSphereNorm,T}

"""
    p = legendre(norm::AbstractLegendreNorm, l::Integer, x::Number)
    P = legendre.(norm::AbstractLegendreNorm, l, x)

Computes the associated Legendre polynomial assuming the order ``m = 0``;
equivalent to `legendre(norm, l, 0, x)` and `legendre.(norm, l, 0, x)`.
"""
@inline function legendre(norm::AbstractLegendreNorm, l::Integer, x::Number)
    return legendre(norm, l, 0, x)
end

"""
    p = Pl(l::Integer, x::Number)

Computes the Legendre polynomials using unit normalization and for degree ``m = 0``;
equivalent to `p = legendre(LegendreUnitNorm(), l, 0, x)`.
"""
@inline Pl(l::Integer, x::Number) = legendre(LegendreUnitNorm(), l, 0, x)

"""
    p = Plm(l::Integer, m::Integer, x::Number)

Computes the associated Legendre polynomials using unit normalization;
equivalent to `p = legendre(LegendreUnitNorm(), l, m, x)`.
"""
@inline Plm(l::Integer, m::Integer, x::Number) = legendre(LegendreUnitNorm(), l, m, x)

"""
    λ = λlm(l::Integer, m::Integer, x::Number)

Computes the associated Legendre polynomials using spherical-harmonic normalization;
equivalent to `λ = legendre(LegendreSphereNorm(), l, m, x)`.
"""
@inline λlm(l::Integer, m::Integer, x::Number) = legendre(LegendreSphereNorm(), l, m, x)

"""
    Pl!(P, l::Integer, x)

Fills the array `P` with the unit-normalized Legendre polynomial values ``P_ℓ(x)`` for fixed
order ``m = 0``;
equivalent to `legendre!(LegendreUnitNorm(), P, l, 0, x)`.
"""
@inline Pl!(P, l::Integer, x) =
    legendre!(LegendreUnitNorm(), P, l, 0, x)

"""
    Plm!(P, l::Integer, m::Integer, x)

Fills the array `P` with the unit-normalized associated Legendre polynomial values
``P_ℓ^m(x)``;
equivalent to `legendre!(LegendreUnitNorm(), P, l, m, x)`.
"""
@inline Plm!(P, l::Integer, m::Integer, x) =
    legendre!(LegendreUnitNorm(), P, l, m, x)


"""
    λlm!(Λ, l::Integer, m::Integer, x)

Fills the array `Λ` with the spherical-harmonic normalized associated Legendre polynomial
values ``λ_ℓ^m(x)``;
equivalent to `legendre!(LegendreSphereNorm(), P, l, m, x)`.
"""
@inline λlm!(Λ, l::Integer, m::Integer, x) =
    legendre!(LegendreSphereNorm(), Λ, l, m, x)

