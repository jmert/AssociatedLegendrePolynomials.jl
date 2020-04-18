# Implements the legendre interface for the unit normalization case

"""
    struct LegendreUnitNorm <: AbstractLegendreNorm end

Trait type denoting the unit normalization of the associated Legendre polynomials.
"""
struct LegendreUnitNorm <: AbstractLegendreNorm end

@inline function
Plm_00(::LegendreUnitNorm, ::Type{T}) where T
    return one(T)
end

@inline function
Plm_μ(::LegendreUnitNorm, ::Type{T}, m::Integer) where T
    return convert(T, 2m - 1)
end

@inline function
Plm_ν(::LegendreUnitNorm, ::Type{T}, m::Integer) where T
    return convert(T, 2m + 1)
end

@inline function
Plm_α(::LegendreUnitNorm, ::Type{T}, l::Integer, m::Integer) where T
    return convert(T, 2l - 1) * inv(convert(T, l - m))
end

@inline function
Plm_β(::LegendreUnitNorm, ::Type{T}, l::Integer, m::Integer) where T
    return convert(T, l + m - 1) * inv(convert(T, l - m))
end
