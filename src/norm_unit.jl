# Implements the legendre interface for the unit normalization case

"""
    struct LegendreUnitNorm <: AbstractLegendreNorm end

Trait type denoting the unit normalization of the associated Legendre polynomials.
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
