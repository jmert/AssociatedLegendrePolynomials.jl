# Implements the legendre interface for the unit normalization case

"""
    struct LegendreSphereNorm <: AbstractLegendreNorm end

Trait type denoting the spherical-harmonic normalization of the associated Legendre
polynomials.
"""
struct LegendreSphereNorm <: AbstractLegendreNorm end

@inline function
Plm_00(::LegendreSphereNorm, ::Type{T}) where T
    # comparing this against
    #   convert(T, inv(sqrt(4*convert(BigFloat, π))))
    # shows that this is exact within Float64 precision
    return inv(sqrt(4 * convert(T, π)))
end

@inline function
Plm_μ(::LegendreSphereNorm, ::Type{T}, l::Integer) where T
    return @fastmath(sqrt)(one(T) + inv(convert(T, 2l)))
end

@inline function
Plm_ν(::LegendreSphereNorm, ::Type{T}, l::Integer) where T
    return @fastmath(sqrt)(convert(T, 2l + 1))
end

@inline function
Plm_α(::LegendreSphereNorm, ::Type{T}, l::Integer, m::Integer) where T
    lT = convert(T, l)
    mT = convert(T, m)
    # Write factors in two pieces to make compiler's job easier. In the case where
    # both Plm_α and Plm_β are called and inlined, the next line from both functions
    # should be merged and shared.
    fac1 = (2lT + 1) / ((2lT - 3) * (lT^2 - mT^2))
    fac2 = 4*(lT - 1)^2 - 1
    @fastmath return sqrt(fac1 * fac2)
end

@inline function
Plm_β(::LegendreSphereNorm, ::Type{T}, l::Integer, m::Integer) where T
    lT = convert(T, l)
    mT = convert(T, m)
    # Write factors in two pieces to make compiler's job easier. In the case where
    # both Plm_α and Plm_β are called and inlined, the next line from both functions
    # should be merged and shared.
    fac1 = (2lT + 1) / ((2lT - 3) * (lT^2 - mT^2))
    fac2 = (lT - 1)^2 - mT^2
    @fastmath return sqrt(fac1 * fac2)
end

# Extra functions

"""
    N = Nlm([T=Float64], l, m)

Computes the normalization constant
```math
    N_ℓ^m ≡ \\sqrt{\\frac{2ℓ+1}{4π} \\frac{(ℓ-m)!}{(ℓ+m)!}}
```
which defines the Spherical Harmonic normalized functions ``λ_ℓ^m(x)`` in
terms of the standard unit normalized ``P_ℓ^m(x)``
```math
    λ_ℓ^m(x) ≡ N_ℓ^m P_ℓ^m(x)
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
