"""
    LegendreUnitCoeff{T}

Precomputed recursion relation coefficients for the standard unit
normalization. Alias for `LegendreNormCoeff{LegendreUnitNorm,T}`.
"""
LegendreUnitCoeff{T} = LegendreNormCoeff{LegendreUnitNorm,T}

"""
    LegendreOrthoCoeff{T}

Table type of precomputed recursion relation coefficients for the orthonormal
normalization. Alias for `LegendreNormCoeff{LegendreOrthoNorm,T}`.
"""
LegendreOrthoCoeff{T} = LegendreNormCoeff{LegendreOrthoNorm,T}

"""
    LegendreSphereCoeff{T}

Table type of precomputed recursion relation coefficients for the spherical
harmonic normalization. Alias for `LegendreNormCoeff{LegendreSphereNorm,T}`.
"""
LegendreSphereCoeff{T} = LegendreNormCoeff{LegendreSphereNorm,T}

"""
    p = Plm(l, m, x)

Computes the associated Legendre polynomials using unit normalization;
equivalent to `p = legendre(LegendreUnitNorm(), l, m, x)`.
"""
const Plm = LegendreUnitNorm()

"""
    Plm!(P, l, m, x)

Fills the array `P` with the unit-normalized associated Legendre polynomial values
``P_ℓ^m(x)``;
equivalent to `legendre!(LegendreUnitNorm(), P, l, m, x)`.
"""
const Plm! = LegendreUnitNorm()

"""
    λ = λlm(l, m, x)

Computes the associated Legendre polynomials using spherical-harmonic normalization;
equivalent to `λ = legendre(LegendreSphereNorm(), l, m, x)`.
"""
const λlm = LegendreSphereNorm()

"""
    λlm!(Λ, l, m, x)

Fills the array `Λ` with the spherical-harmonic normalized associated Legendre polynomial
values ``λ_ℓ^m(x)``;
equivalent to `legendre!(LegendreSphereNorm(), P, l, m, x)`.
"""
const λlm! = LegendreSphereNorm()
