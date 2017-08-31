"""
Collections of functions which compute the associated Legendre functions.

Based on implementation described in Limpanuparb and Milthorpe (2014)
*“Associated Legendre Polynomials and Spherical Harmonics Computation for
Chemistry Applications”* arXiv:1410.1748v1
"""
module Legendre

# Normalization trait and cache types
export
    AbstractLegendreNorm, LegendreNormCoeff,
    LegendreUnitNorm,     LegendreUnitCoeff,
    LegendreSphereNorm,   LegendreSphereCoeff

# Generic computation functions
export LegendreP, LegendreP!

# Specific computation functions
export Pl, Pl!, Plm, Plm!, Nlm, λlm, λlm!

import Base.@boundscheck, Base.@propagate_inbounds



"""
    abstract type AbstractLegendreNorm end

Abstract supertype for normalization conditions of the Associated Legendre polynomials.
"""
abstract type AbstractLegendreNorm end

"""
    struct LegendreUnitNorm <: AbstractLegendreNorm end
"""
struct LegendreUnitNorm <: AbstractLegendreNorm end

"""
    struct LegendreSphereNorm <: AbstractLegendreNorm end
"""
struct LegendreSphereNorm <: AbstractLegendreNorm end

"""
    type LegendreNormCoeff{N<:AbstractLegendreNorm,T<:Real} <: AbstractLegendreNorm

Precomputed recursion relation coefficients for the normalization `N` and value type
`T`.

# Example
```jldoctest
julia> LegendreNormCoeff{LegendreSphereNorm,Float64}(1)
CMB.Legendre.LegendreNormCoeff{CMB.Legendre.LegendreSphereNorm,Float64} for lmax = 1 with coefficients:
    μ: [0.0, 1.22474]
    ν: [0.0, 1.73205]
    α: [0.0 1.93649; 1.73205 0.0]
    β: [0.0 1.11803; -0.0 0.0]
```
"""
struct LegendreNormCoeff{N<:AbstractLegendreNorm,T<:Real} <: AbstractLegendreNorm
    μ::Vector{T}
    ν::Vector{T}
    α::Matrix{T}
    β::Matrix{T}

    function LegendreNormCoeff{N,T}(lmax::Integer) where {N,T}
        lmax ≥ 0 || throw(DomainError())

        μ = zeros(T, lmax+1)
        ν = zeros(T, lmax+1)
        α = zeros(T, lmax+1, lmax+1)
        β = zeros(T, lmax+1, lmax+1)

        @inbounds for l in 0:lmax
            μ[l+2] = Plm_μ(N(), T, l+1)
            ν[l+2] = Plm_ν(N(), T, l+1)

            for m in 0:l
                α[l+2,m+1] = Plm_α(N(), T, l+1, m)
                β[l+2,m+1] = Plm_β(N(), T, l+1, m)
            end
        end

        return new(μ, ν, α, β)
    end
end

"""
    LegendreUnitCoeff{T}

Precomputed recursion relation coefficients for the standard unit
normalization. Alias for `LegendreNormCoeff{LegendreUnitNorm,T}`.
"""
LegendreUnitCoeff = LegendreNormCoeff{LegendreUnitNorm}

"""
    LegendreSphereCoeff{T}

Table type of precomputed recursion relation coefficients for the spherical
harmonic normalization. Alias for `LegendreNormCoeff{LegendreSphereNorm,T}`.
"""
LegendreSphereCoeff = LegendreNormCoeff{LegendreSphereNorm}

# Improve printing somewhat
Base.show(io::IO, norm::LegendreNormCoeff{N,T}) where {N,T} =
    print(io, LegendreNormCoeff, "{$N,$T}")
function Base.show(io::IO, ::MIME"text/plain", P::LegendreNormCoeff)
    println(io, P, " for lmax = $(length(P.μ)-1) with coefficients:")
    println(io, "    μ: ", P.μ)
    println(io, "    ν: ", P.ν)
    println(io, "    α: ", P.α)
    println(io, "    β: ", P.β)
end

"""
    Plm_00(::N, ::Type{T}) where {N<:AbstractLegendreNorm, T<:Real}

Returns the initial condition ``P_0^0(x)`` for the associated Legendre recursions based
on the normalization choice `N` for numeric type `T`.
"""
function Plm_00 end

"""
    Plm_μ(norm::N, ::Type{T}, l::Integer) where {N<:AbstractLegendreNorm, T<:Real}

Returns the coefficient ``μ_ℓ`` for the single-term recursion relation
```math
    P_{ℓ+1}^{ℓ+1}(x) = -μ_{ℓ+1} \\sqrt{1-x^2} P_ℓ^ℓ(x)
```
where ``μ_ℓ`` is appropriate for the choice of normalization `N`.
"""
function Plm_μ end

"""
    Plm_ν(norm::N, ::Type{T}, l::Integer) where {N<:AbstractLegendreNorm, T<:Real}

Returns the coefficient ``ν_ℓ`` for the single-term recursion relation
```math
    P_{ℓ+1}^ℓ(x) = ν_{ℓ+1} x P_ℓ^ℓ(x)
```
where ``ν_ℓ`` is appropriate for the choice of normalization `N`.
"""
function Plm_ν end

"""
    Plm_α(norm::N, ::Type{T}, l::Integer, m::Integer) where {N<:AbstractLegendreNorm, T<:Real}

Returns the coefficient ``α_ℓ^m`` for the two-term recursion relation
```math
    P_{ℓ+1}^{m}(x) = α_{ℓ+1}^m x P_ℓ^m(x) - β_{ℓ+1}^m P_{ℓ-1}^m(x)
```
where ``α_ℓ^m`` is appropriate for the choice of normalization `N`.
"""
function Plm_α end

"""
    Plm_β(norm::N, ::Type{T}, l::Integer, m::Integer) where {N<:AbstractLegendreNorm, T<:Real}

Returns the coefficient ``β_ℓ^m`` for the two-term recursion relation
```math
    P_{ℓ+1}^{m}(x) = α_{ℓ+1}^m x P_ℓ^m(x) - β_{ℓ+1}^m P_{ℓ-1}^m(x)
```
where ``β_ℓ^m`` is appropriate for the choice of normalization `N`.
"""
function Plm_β end

@inline function
Plm_00(::LegendreUnitNorm, ::Type{T}) where {T<:Real}
    return one(T)
end

@inline function
Plm_00(::LegendreSphereNorm, ::Type{T}) where {T<:Real}
    # comparing this against
    #   convert(T, inv(sqrt(4*convert(BigFloat, π))))
    # shows that this is exact within Float64 precision
    return convert(T, inv(sqrt(4π)))
end

@inline function
Plm_00(::LegendreNormCoeff{N,T}, ::Type{T}
                       ) where {N<:AbstractLegendreNorm, T<:Real}
    return Plm_00(N(), T)
end

@inline function
Plm_μ(::LegendreUnitNorm, ::Type{T}, l::Integer) where {T<:Real}
    return convert(T, 2l - 1)
end

@inline function
Plm_μ(::LegendreSphereNorm, ::Type{T}, l::Integer) where {T<:Real}
    return sqrt(one(T) + inv(convert(T, 2l)))
end

@propagate_inbounds function
Plm_μ(norm::LegendreNormCoeff, ::Type{T}, l::Integer) where {T<:Real}
    return norm.μ[l+1]
end

@inline function
Plm_ν(::LegendreUnitNorm, ::Type{T}, l::Integer) where {T<:Real}
    return convert(T, 2l - 1)
end

@inline function
Plm_ν(::LegendreSphereNorm, ::Type{T}, l::Integer) where {T<:Real}
    return sqrt(convert(T, 2l + 1))
end

@propagate_inbounds function
Plm_ν(norm::LegendreNormCoeff, ::Type{T}, l::Integer) where {T<:Real}
    return norm.ν[l+1]
end

@inline function
Plm_α(::LegendreUnitNorm, ::Type{T}, l::Integer, m::Integer) where T<:Real
    return convert(T, 2l-1) * inv(convert(T, l-m))
end

@inline function
Plm_α(::LegendreSphereNorm, ::Type{T}, l::Integer, m::Integer) where T<:Real
    lT = convert(T, l)
    mT = convert(T, m)
    # Write factors in two pieces to make compiler's job easier. In the case where
    # both Plm_α and Plm_β are called and inlined, the next line from both functions
    # should be merged and shared.
    fac1 = (2lT + 1) * inv((2lT - 3) * (lT*lT - mT*mT))
    # Explicitly perform ((lT-1)*(lT-1)) as a group to allow compiler to reuse (ℓ-1)² if
    # Plm_β is being computed, too. (Remember, floating point math is non-commutative,
    # so the compiler isn't free to change (4*(lT-1))*(lT-1) to 4*((lT-1)*(lT-1)).)
    fac2 = 4*((lT-1)*(lT-1)) - 1
    @fastmath return sqrt(fac1 * fac2)
end

@propagate_inbounds function
Plm_α(norm::LegendreNormCoeff, ::Type{T}, l::Integer, m::Integer
        ) where {T<:Real}
    return norm.α[l+1,m+1]
end

@inline function
Plm_β(::LegendreUnitNorm, ::Type{T}, l::Integer, m::Integer
        ) where T<:Real
    return convert(T, l+m-1) * inv(convert(T, l-m))
end

@inline function
Plm_β(::LegendreSphereNorm, ::Type{T}, l::Integer, m::Integer) where T<:Real
    lT = convert(T, l)
    mT = convert(T, m)
    # Write factors in two pieces to make compiler's job easier. In the case where
    # both Plm_α and Plm_β are called and inlined, the next line from both functions
    # should be merged and shared.
    fac1 = (2lT + 1) * inv((2lT - 3) * (lT*lT - mT*mT))
    fac2 = (lT-1)*(lT-1) - mT*mT
    @fastmath return sqrt(fac1 * fac2)
end

@propagate_inbounds function
Plm_β(norm::LegendreNormCoeff, ::Type{T}, l::Integer, m::Integer) where {T<:Real}
    return norm.β[l+1,m+1]
end

# Named recursion relations

@propagate_inbounds function
_Plm_1term_raise_lm(norm::N, l::Integer, x::T, y::T, plm::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    μ = Plm_μ(norm, T, l+1)
    return -μ * y * plm
end
@propagate_inbounds function
_Plm_1term_raise_l(norm::N, l::Integer, x::T, plm1m::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    ν = Plm_ν(norm, T, l+1)
    return ν * x * plm1m
end
@propagate_inbounds function
_Plm_2term_raise_l(norm::N, l::Integer, m::Integer, x::T, plm1m::T, plm2m::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    α = Plm_α(norm, T, l+1, m)
    β = Plm_β(norm, T, l+1, m)
    return α * x * plm1m - β * plm2m
end

# Bounds checks, both mathematical domain and array access (if applicable)

@noinline function
_chkbounds_Pl(norm::N, l::Integer, x::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    (0 ≤ l) || throw(DomainError())
end
@noinline function
_chkbounds_Pl(norm::LegendreNormCoeff{N,T}, l::Integer, x::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    (0 ≤ l) || throw(DomainError())
    (l ≤ length(norm.μ)) || throw(BoundsError())
end
@noinline function
_chkbounds_Plm(norm::N, l::Integer, m::Integer, x::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    (0 ≤ l && 0 ≤ m ≤ l) || throw(DomainError())
end
@noinline function
_chkbounds_Plm(norm::LegendreNormCoeff{N,T}, l::Integer, m::Integer, x::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    (0 ≤ l && 0 ≤ m ≤ l) || throw(DomainError())
    (l ≤ length(norm.μ)) || throw(BoundsError())
end

@noinline function
_chkbounds_Pl!(norm::N, Λ::AbstractVector{T}, lmax::Integer, m::Integer, x::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    (0 ≤ lmax && 0 ≤ m ≤ lmax) || throw(DomainError())
    (size(Λ,1)≥lmax+1) || throw(DimensionMismatch())
end
@noinline function
_chkbounds_Pl!(norm::LegendreNormCoeff{N,T}, Λ::AbstractVector{T}, lmax::Integer,
        m::Integer, x::T) where {N<:AbstractLegendreNorm, T<:Real}
    (0 ≤ lmax && 0 ≤ m ≤ lmax) || throw(DomainError())
    (lmax ≤ length(norm.μ)) || throw(BoundsError())
    (size(Λ,1)≥lmax+1) || throw(DimensionMismatch())
end

@noinline function
_chkbounds_Plm!(norm::N, Λ::AbstractMatrix{T}, lmax::Integer, x::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    (0 ≤ lmax) || throw(DomainError())
    (size(Λ,1)≥lmax+1 && size(Λ,2)≥lmax+1) || throw(DimensionMismatch())
end
@noinline function
_chkbounds_Plm!(norm::LegendreNormCoeff{N,T}, Λ::AbstractMatrix{T}, lmax::Integer,
        x::T) where {N<:AbstractLegendreNorm, T<:Real}
    (0 ≤ lmax) || throw(DomainError())
    (lmax ≤ length(norm.μ)) || throw(BoundsError())
    (size(Λ,1)≥lmax+1 && size(Λ,2)≥lmax+1) || throw(DimensionMismatch())
end

"""
Computes the Associated Legendre polynomials ``P_ℓ^m(x)``.

    LegendreP(norm, l, x) -> P
    LegendreP(norm, l, m, x) -> P

    LegendreP!(norm, P::Vector, lmax, m, x) -> P
    LegendreP!(norm, P::Matrix, lmax, x) -> P
"""
function LegendreP end, function LegendreP! end

function LegendreP(norm::N, l::Integer, x::T
                  ) where {N<:AbstractLegendreNorm, T<:Real}
    @boundscheck _chkbounds_Pl(norm, l, x)
    @inbounds begin
        pl = Plm_00(norm, T)
        l == 0 && return pl

        # now boost one in l to P_1 using a single-term recurrence
        plp1 = _Plm_1term_raise_l(norm, 0, x, pl)
        plm1 = pl
        pl = plp1

        # then finish by iterating to P_l^m using the two-term recurrence
        for n in 1:(l-1)
            plp1 = _Plm_2term_raise_l(norm, n, 0, x, pl, plm1)
            plm1 = pl
            pl = plp1
        end

        return pl
    end
end

function LegendreP!(norm::N, P::AbstractVector{T}, lmax::Integer, m::Integer, x::T
                   ) where {N<:AbstractLegendreNorm, T<:Real}
    @boundscheck _chkbounds_Pl!(norm, P, lmax, m, x)
    @inbounds begin
        pl = Plm_00(norm, T)
        P[1] = pl
        lmax == 0 && return P

        # Iterate along the main diagonal until we reach the target m
        y = sqrt(one(T) - x*x)
        for n in 0:(m-1)
            plp1 = _Plm_1term_raise_lm(norm, n, x, y, pl)
            pl = plp1
            P[n+1] = zero(T)
        end
        P[m+1] = pl

        # First step is to boost one in l to P_{m+1}^m using a single-term
        # recurrence
        plp1 = _Plm_1term_raise_l(norm, m, x, pl)
        P[m+2] = plp1
        plm1 = pl
        pl = plp1

        # Then finish by iterating to P_lmax^m using the two-term recurrence
        for l in (m+1):(lmax-1)
            plp1 = _Plm_2term_raise_l(norm, l, m, x, pl, plm1)
            P[l+2] = plp1
            plm1 = pl
            pl = plp1
        end # for l

        return P
    end
end

function LegendreP(norm::N, l::Integer, m::Integer, x::T
                  ) where {N<:AbstractLegendreNorm, T<:Real}
    @boundscheck _chkbounds_Plm(norm, l, m, x)
    @inbounds begin
        pl = Plm_00(norm, T)
        l == 0 && return pl

        # iterate along the diagonal to P_m^m
        y = sqrt(one(T) - x*x)
        for n in 0:(m-1)
            plp1 = _Plm_1term_raise_lm(norm, n, x, y, pl)
            pl = plp1
        end
        # return if an on-diagonal term was requested
        l == m && return pl

        # now boost one in l to P_{m+1}^m using a single-term recurrence
        plp1 = _Plm_1term_raise_l(norm, m, x, pl)
        plm1 = pl
        pl = plp1

        # then finish by iterating to P_l^m using the two-term recurrence
        for n in (m+1):(l-1)
            plp1 = _Plm_2term_raise_l(norm, n, m, x, pl, plm1)
            plm1 = pl
            pl = plp1
        end

        return pl
    end
end

function LegendreP!(norm::N, Λ::AbstractMatrix{T}, lmax::Integer, x::T
                   ) where {N<:AbstractLegendreNorm, T<:Real}
    @boundscheck _chkbounds_Plm!(norm, Λ, lmax, x)
    @inbounds begin
        Λ[1,1] = Plm_00(norm, T)
        lmax == 0 && return Λ

        # Fill in the main diagonal first
        y = sqrt(one(T) - x*x)
        pl = Λ[1,1]
        for l in 0:(lmax-1)
            plp1 = _Plm_1term_raise_lm(norm, l, x, y, pl)
            Λ[l+2,l+2] = plp1
            pl = plp1
        end

        # Outer loop runs over the m's
        for m in 0:lmax
            # First step is to boost one in l to P_{m+1}^m using a single-term
            # recurrence
            pl = Λ[m+1,m+1]
            plp1 = _Plm_1term_raise_l(norm, m, x, pl)
            Λ[m+2,m+1] = plp1
            plm1 = pl
            pl = plp1

            # Then finish by iterating to P_lmax^m using the two-term recurrence
            for l in (m+1):(lmax-1)
                plp1 = _Plm_2term_raise_l(norm, l, m, x, pl, plm1)
                Λ[l+2,m+1] = plp1
                plm1 = pl
                pl = plp1
            end # for l
        end # for m

        return Λ
    end
end

"""
Computes the Legendre polynomials ``P_ℓ(x)``.

    Pl(l, x) -> P
    Pl!(P::Vector, lmax, x) -> P

"""
function Pl end, function Pl! end

"""
Computes the Associated Legendre polynomial ``P_ℓ^m(x)`` of order ``ℓ`` and degree
``m``.

    Plm(l, m, x) -> P
    Plm!(P::Matrix, lmax, x) -> P
"""
function Plm end, function Plm! end

"""
Computes the spherical harmonic pre-normalized Associated Legendre polynomial
``λ_ℓ^m(x)`` of order ``ℓ`` and degree ``m``.

    λlm(l, m, x) -> Λ
    λlm!(Λ::Matrix, lmax, x) -> Λ
"""
function λlm end, function λlm! end

@inline Pl(l::Integer, x::T) where {T<:Real} = LegendreP(LegendreUnitNorm(), l, x)
@inline Pl!(P::AbstractVector{T}, lmax::Integer, x::T) where {T<:Real} =
    LegendreP!(LegendreUnitNorm(), P, lmax, 0, x)

@inline Plm(l::Integer, m::Integer, x::T) where {T<:Real} =
    LegendreP(LegendreUnitNorm(),   l, m, x)
@inline λlm(l::Integer, m::Integer, x::T) where {T<:Real} =
    LegendreP(LegendreSphereNorm(), l, m, x)

@inline Plm!(P::AbstractMatrix{T}, lmax::Integer, x::T) where {T<:Real} =
    LegendreP!(LegendreUnitNorm(), P, lmax, x)
@inline λlm!(Λ::AbstractMatrix{T}, lmax::Integer, x::T) where {T<:Real} =
    LegendreP!(LegendreSphereNorm(), Λ, lmax, x)

"""
    Nlm(::Type{<:Real}, l, m) -> N

Computes the normalization constant
```math
    N_ℓ^m ≡ \\sqrt{\\frac{2ℓ+1}{4π} \\frac{(ℓ-m)!}{(ℓ+m)!}}
```
which defines the Spherical Harmonic normalized functions ``λ_ℓ^m(x)`` in
terms of the standard normalized ``P_ℓ^m(x)``:
```math
    λ_ℓ^m(x) ≡ N_ℓ^m P_ℓ^m(x)
```
See also [`Plm`](@ref) and [`λlm`](@ref).
"""
function Nlm end

function Nlm(::Type{T}, l::Integer, m::Integer) where T<:Real
    fac1 = one(T)
    for ii in (l-m+1):(l+m)
        fac1 *= convert(T, ii)
    end
    num = convert(T, 2l+1)
    den = convert(T, 4π)
    return sqrt( num * inv(den) * inv(fac1) )
end

end # module Legendre
