# Implements the legendre interface for the unit normalization case

"""
    struct LegendreUnitNorm <: AbstractLegendreNorm end

Trait type denoting the unnormalized associated Legendre functions ``P_\\ell^m(x)``
which solve the colatitude ``\\theta`` part of Laplace's equation in spherical coordinates
where ``x=\\cos(\\theta)``.
The degree ``\\ell=0``, order ``m=0`` constant function is normalized to be unity:
```math
    P_0^0(x) = 1
```
"""
struct LegendreUnitNorm <: AbstractLegendreNorm end

@inline function
initcond(::LegendreUnitNorm, ::Type{T}) where T
    return one(T)
end

@inline function
coeff_μ(::LegendreUnitNorm, ::Type{T}, l::Integer) where T
    return convert(T, 2l - 1)
end

@inline function
coeff_ν(::LegendreUnitNorm, ::Type{T}, l::Integer) where T
    return convert(T, 2l - 1)
end

@inline function
coeff_α(::LegendreUnitNorm, ::Type{T}, l::Integer, m::Integer) where T
    return convert(T, 2l - 1) * inv(convert(T, l - m))
end

@inline function
coeff_β(::LegendreUnitNorm, ::Type{T}, l::Integer, m::Integer) where T
    return convert(T, l + m - 1) * inv(convert(T, l - m))
end
