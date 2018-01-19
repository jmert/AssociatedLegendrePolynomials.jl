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
export legendre, legendre!

# Specific computation functions
export Pl, Pl!, Plm, Plm!, Nlm, λlm, λlm!

import Base.@boundscheck, Base.@propagate_inbounds



"""
    abstract type AbstractLegendreNorm end

Abstract supertype for normalization conditions of the Associated Legendre polynomials.

# Example
```jldoctest
julia> subtypes(AbstractLegendreNorm)
3-element Array{Union{DataType, UnionAll},1}:
 CMB.Legendre.LegendreNormCoeff
 CMB.Legendre.LegendreSphereNorm
 CMB.Legendre.LegendreUnitNorm
```
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
    struct LegendreNormCoeff{N<:AbstractLegendreNorm,T<:Real} <: AbstractLegendreNorm

Precomputed recursion relation coefficients for the normalization `N` and value type
`T`.

# Example
```jldoctest
julia> LegendreNormCoeff{LegendreSphereNorm,Float64}(1)
CMB.Legendre.LegendreNormCoeff{CMB.Legendre.LegendreSphereNorm,Float64} for lmax = 1, mmax = 1 with coefficients:
    μ: [0.0, 1.22474]
    ν: [1.73205, 2.23607]
    α: [0.0 0.0; 1.73205 0.0]
    β: [0.0 0.0; -0.0 0.0]
```
"""
struct LegendreNormCoeff{N<:AbstractLegendreNorm,T<:Real} <: AbstractLegendreNorm
    μ::Vector{T}
    ν::Vector{T}
    α::Matrix{T}
    β::Matrix{T}

    function LegendreNormCoeff{N,T}(lmax::Integer, mmax::Integer) where {N,T}
        (lmax ≥ 0 && 0 ≤ mmax ≤ lmax) || throw(DomainError())

        μ = zeros(T, mmax+1)
        ν = zeros(T, mmax+1)
        α = zeros(T, lmax+1, mmax+1)
        β = zeros(T, lmax+1, mmax+1)

        @inbounds for m in 0:mmax
            μ[m+1] = m == 0 ? zero(T) : Plm_μ(N(), T, m)
            ν[m+1] = Plm_ν(N(), T, m)

            for l in (m+1):lmax
                α[l+1,m+1] = Plm_α(N(), T, l, m)
                β[l+1,m+1] = Plm_β(N(), T, l, m)
            end
        end

        return new(μ, ν, α, β)
    end
end

LegendreNormCoeff{N,T}(lmax::Integer) where {N,T} = LegendreNormCoeff{N,T}(lmax, lmax)

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
function Base.show(io::IO, ::MIME"text/plain", N::LegendreNormCoeff)
    lmax,mmax = size(N.α) .- 1
    println(io, N, " for lmax = $lmax, mmax = $mmax with coefficients:")
    println(io, "    μ: ", N.μ)
    println(io, "    ν: ", N.ν)
    println(io, "    α: ", N.α)
    println(io, "    β: ", N.β)
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
Plm_μ(::LegendreUnitNorm, ::Type{T}, m::Integer) where {T<:Real}
    return convert(T, 2m - 1)
end

@inline function
Plm_μ(::LegendreSphereNorm, ::Type{T}, m::Integer) where {T<:Real}
    return sqrt(one(T) + inv(convert(T, 2m)))
end

@propagate_inbounds function
Plm_μ(norm::LegendreNormCoeff, ::Type{T}, m::Integer) where {T<:Real}
    return norm.μ[m+1]
end

@inline function
Plm_ν(::LegendreUnitNorm, ::Type{T}, m::Integer) where {T<:Real}
    return convert(T, 2m + 1)
end

@inline function
Plm_ν(::LegendreSphereNorm, ::Type{T}, m::Integer) where {T<:Real}
    return sqrt(convert(T, 2m + 3))
end

@propagate_inbounds function
Plm_ν(norm::LegendreNormCoeff, ::Type{T}, m::Integer) where {T<:Real}
    return norm.ν[m+1]
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
_1term_raise_lm(norm::N, m::Integer, x::T, y::T, plm::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    μ = Plm_μ(norm, T, m+1)
    return -μ * y * plm
end
@propagate_inbounds function
_1term_raise_l(norm::N, m::Integer, x::T, plm1m::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    ν = Plm_ν(norm, T, m)
    return ν * x * plm1m
end
@propagate_inbounds function
_2term_raise_l(norm::N, l::Integer, m::Integer, x::T, plm1m::T, plm2m::T
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
    (l < size(norm.α,1)) || throw(BoundsError())
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
    (l < size(norm.α,1) && m < size(norm.α,2)) || throw(BoundsError())
end

@noinline function
_chkbounds_Pl!(norm::N, Λ::AbstractVector{T}, lmax::Integer, m::Integer, x::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    (0 ≤ lmax && 0 ≤ m ≤ lmax) || throw(DomainError())
    (size(Λ,1) ≥ lmax+1) || throw(DimensionMismatch())
end
@noinline function
_chkbounds_Pl!(norm::LegendreNormCoeff{N,T}, Λ::AbstractVector{T}, lmax::Integer,
        m::Integer, x::T) where {N<:AbstractLegendreNorm, T<:Real}
    (0 ≤ lmax && 0 ≤ m ≤ lmax) || throw(DomainError())
    (lmax < size(norm.α,1) && m < size(norm.α,2)) || throw(BoundsError())
    (size(Λ,1) ≥ lmax+1) || throw(DimensionMismatch())
end

@noinline function
_chkbounds_Plm!(norm::N, Λ::AbstractMatrix{T}, lmax::Integer, mmax::Integer, x::T
        ) where {N<:AbstractLegendreNorm, T<:Real}
    (0 ≤ lmax && 0 ≤ mmax ≤ lmax) || throw(DomainError())
    (size(Λ,1) ≥ lmax+1 && size(Λ,2) ≥ mmax+1) || throw(DimensionMismatch())
end
@noinline function
_chkbounds_Plm!(norm::LegendreNormCoeff{N,T}, Λ::AbstractMatrix{T}, lmax::Integer,
        mmax::Integer, x::T) where {N<:AbstractLegendreNorm, T<:Real}
    (0 ≤ lmax && 0 ≤ mmax ≤ lmax) || throw(DomainError())
    (lmax < size(norm.α,1) && mmax < size(norm.α,2)) || throw(BoundsError())
    (size(Λ,1) ≥ lmax+1 && size(Λ,2) ≥ mmax+1) || throw(DimensionMismatch())
end

"""
    p = legendre(norm::AbstractLegendreNorm, l::Integer, x::Real)

Computes the scalar value ``p = N_ℓ P_ℓ(x)``, where ``P_ℓ(x)`` is the Legendre
polynomial of degree `l` at `x` and ``N_ℓ`` is the normalization scheme `norm`.
"""
function legendre(norm::N, l::Integer, x::T
                  ) where {N<:AbstractLegendreNorm, T<:Real}
    @boundscheck _chkbounds_Pl(norm, l, x)
    @inbounds begin
        pl = Plm_00(norm, T)
        l == 0 && return pl

        # now boost one in l to P_1 using a single-term recurrence
        plp1 = _1term_raise_l(norm, 0, x, pl)
        plm1 = pl
        pl = plp1

        # then finish by iterating to P_l^m using the two-term recurrence
        for n in 1:(l-1)
            plp1 = _2term_raise_l(norm, n, 0, x, pl, plm1)
            plm1 = pl
            pl = plp1
        end

        return pl
    end
end

"""
    legendre!(norm::AbstractLegendreNorm, Λ::AbstractVector, lmax::Integer, m::Integer,
              x::Real)

Fills the vector `Λ` with the pre-normalized Legendre polynomial values ``N_ℓ^m P_ℓ^m(x)``
for all degrees `0 ≤ ℓ ≤ lmax` and constant order `m` at `x`, where ``N_ℓ^m`` is the
normalization scheme `norm`.
"""
function legendre!(norm::N, Λ::AbstractVector{T}, lmax::Integer, m::Integer, x::T
                   ) where {N<:AbstractLegendreNorm, T<:Real}
    @boundscheck _chkbounds_Pl!(norm, Λ, lmax, m, x)
    @inbounds begin
        pl = Plm_00(norm, T)
        Λ[1] = pl
        lmax == 0 && return Λ

        # Iterate along the main diagonal until we reach the target m
        y = sqrt(one(T) - x*x)
        for n in 0:(m-1)
            plp1 = _1term_raise_lm(norm, n, x, y, pl)
            pl = plp1
            Λ[n+1] = zero(T)
        end
        Λ[m+1] = pl

        # First step is to boost one in l to P_{m+1}^m using a single-term
        # recurrence
        plp1 = _1term_raise_l(norm, m, x, pl)
        Λ[m+2] = plp1
        plm1 = pl
        pl = plp1

        # Then finish by iterating to P_lmax^m using the two-term recurrence
        for l in (m+1):(lmax-1)
            plp1 = _2term_raise_l(norm, l, m, x, pl, plm1)
            Λ[l+2] = plp1
            plm1 = pl
            pl = plp1
        end # for l

        return Λ
    end
end

"""
    p = legendre(norm::AbstractLegendreNorm, l::Integer, m::Integer, x::Real)

Computes the scalar value ``p = N_ℓ^m P_ℓ^m(x)``, where ``P_ℓ^m(x)`` is the associated
Legendre polynomial of degree `l` and order `m` at `x` and ``N_ℓ^m`` the normalization
scheme `norm`.
"""
function legendre(norm::N, l::Integer, m::Integer, x::T
                  ) where {N<:AbstractLegendreNorm, T<:Real}
    @boundscheck _chkbounds_Plm(norm, l, m, x)
    @inbounds begin
        pl = Plm_00(norm, T)
        l == 0 && return pl

        # iterate along the diagonal to P_m^m
        y = sqrt(one(T) - x*x)
        for n in 0:(m-1)
            plp1 = _1term_raise_lm(norm, n, x, y, pl)
            pl = plp1
        end
        # return if an on-diagonal term was requested
        l == m && return pl

        # now boost one in l to P_{m+1}^m using a single-term recurrence
        plp1 = _1term_raise_l(norm, m, x, pl)
        plm1 = pl
        pl = plp1

        # then finish by iterating to P_l^m using the two-term recurrence
        for n in (m+1):(l-1)
            plp1 = _2term_raise_l(norm, n, m, x, pl, plm1)
            plm1 = pl
            pl = plp1
        end

        return pl
    end
end

"""
    legendre!(norm::AbstractLegendreNorm, P::AbstractMatrix, lmax::Integer, mmax::Integer, x::Real)

Fills the matrix `Λ` with the pre-normalized Legendre polynomial values ``N_ℓ^m P_ℓ^m(x)``
for all degrees `0 ≤ ℓ ≤ lmax` and all orders `0 ≤ m ≤ ℓ` at `x`, where ``N_ℓ^m`` is the
normalization scheme `norm`.
"""
function legendre!(norm::N, Λ::AbstractMatrix{T}, lmax::Integer, mmax::Integer, x::T
                   ) where {N<:AbstractLegendreNorm, T<:Real}
    @boundscheck _chkbounds_Plm!(norm, Λ, lmax, mmax, x)
    @inbounds begin
        Λ[1,1] = Plm_00(norm, T)
        lmax == 0 && return Λ

        # Fill in the main diagonal first
        y = sqrt(one(T) - x*x)
        pl = Λ[1,1]
        for m in 0:(mmax-1)
            plp1 = _1term_raise_lm(norm, m, x, y, pl)
            Λ[m+2,m+2] = plp1
            pl = plp1
        end

        # Outer loop runs over the m's
        for m in 0:min(mmax, lmax-1)
            # First step is to boost one in l to P_{m+1}^m using a single-term
            # recurrence
            pl = Λ[m+1,m+1]
            plp1 = _1term_raise_l(norm, m, x, pl)
            Λ[m+2,m+1] = plp1
            plm1 = pl
            pl = plp1

            # Then finish by iterating to P_lmax^m using the two-term recurrence
            for l in (m+1):(lmax-1)
                plp1 = _2term_raise_l(norm, l, m, x, pl, plm1)
                Λ[l+2,m+1] = plp1
                plm1 = pl
                pl = plp1
            end # for l
        end # for m

        return Λ
    end
end

# Make coefficient cache objects callable with similar syntax as the legendre[!] functions
@propagate_inbounds function (norm::LegendreNormCoeff)(l::Integer, x::Real)
    return legendre(norm, l, x)
end
@propagate_inbounds function (norm::LegendreNormCoeff)(l::Integer, m::Integer, x::Real)
    return legendre(norm, l, m, x)
end
@propagate_inbounds function (norm::LegendreNormCoeff)(Λ::AbstractVector{T}, m::Integer, x::T) where T
    lmax = size(norm.α,1) - 1
    return legendre!(norm, Λ, lmax, m, x)
end
@propagate_inbounds function (norm::LegendreNormCoeff)(Λ::AbstractMatrix{T}, x::T) where T
    lmax,mmax = size(norm.α) .- 1
    return legendre!(norm, Λ, lmax, mmax, x)
end


"""
    p = Pl(l::Integer, x::Real)

Computes the scalar value ``p = P_ℓ(x)``, where ``P_ℓ(x)`` is the Legendre polynomial of
degree `l` at `x`.
"""
@inline Pl(l::Integer, x::T) where {T<:Real} = legendre(LegendreUnitNorm(), l, x)

"""
    p = Plm(l::Integer, m::Integer, x::Real)

Computes the scalar value ``p = P_ℓ^m(x)``, where ``P_ℓ^m(x)`` is the associated Legendre
polynomial of degree `l` and order `m` at `x`.
"""
@inline Plm(l::Integer, m::Integer, x::T) where {T<:Real} =
    legendre(LegendreUnitNorm(), l, m, x)

"""
    λ = λlm(l::Integer, m::Integer, x::Real)

Computes the scalar value ``λ = λ_ℓ^m(x)``, where ``λ_ℓ^m(x)`` is the spherical-harmonic
normalized associated Legendre polynomial of degree `l` and order `m` at `x`.
"""
@inline λlm(l::Integer, m::Integer, x::T) where {T<:Real} =
    legendre(LegendreSphereNorm(), l, m, x)

"""
    Pl!(P::AbstractVector, lmax::Integer, x::Real)

Fills the vector `P` with the Legendre polynomial values ``P_ℓ(x)`` for all degrees
`0 ≤ ℓ ≤ lmax` at `x`.
"""
@inline Pl!(P::AbstractVector{T}, lmax::Integer, x::T) where {T<:Real} =
    legendre!(LegendreUnitNorm(), P, lmax, 0, x)

"""
    Plm!(P::AbstractVector, lmax::Integer, m::Integer, x::Real)

Fills the vector `P` with the Legendre polynomial values ``P_ℓ^m(x)`` for all degrees
`0 ≤ ℓ ≤ lmax` and constant order `m` at `x`.
"""
@inline Plm!(P::AbstractVector{T}, lmax::Integer, m::Integer, x::T) where {T<:Real} =
    legendre!(LegendreUnitNorm(), P, lmax, m, x)


"""
    Plm!(P::AbstractMatrix, lmax::Integer, mmax::Integer, x::Real)

Fills the lower triangle of the matrix `P` with the associated Legendre polynomial values
``P_ℓ^m(x)`` for all degrees `0 ≤ ℓ ≤ lmax` and all orders `0 ≤ m ≤ ℓ` at `x`.
"""
@inline Plm!(P::AbstractMatrix{T}, lmax::Integer, mmax::Integer, x::T) where {T<:Real} =
    legendre!(LegendreUnitNorm(), P, lmax, mmax, x)

"""
    λlm!(Λ::AbstractVector, lmax::Integer, m::Integer, x::Real)

Fills the vector `Λ` with the spherical harmonic normalized associated Legendre polynomial
values ``λ_ℓ^m(x)`` for all degrees `0 ≤ ℓ ≤ lmax` and constant order `m` at `x`.
"""
@inline λlm!(Λ::AbstractVector{T}, lmax::Integer, m::Integer, x::T) where {T<:Real} =
    legendre!(LegendreSphereNorm(), Λ, lmax, m, x)

"""
    λlm!(Λ::AbstractMatrix, lmax::Integer, mmax::Integer, x::Real)

Fills the lower triangle of the matrix `Λ` with the spherical harmonic normalized associated
Legendre polynomial values ``Λ_ℓ^m(x)`` for all degrees `0 ≤ ℓ ≤ lmax` and all orders
`0 ≤ m ≤ ℓ` at `x`.
"""
@inline λlm!(Λ::AbstractMatrix{T}, lmax::Integer, mmax::Integer, x::T) where {T<:Real} =
    legendre!(LegendreSphereNorm(), Λ, lmax, mmax, x)

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
function Nlm(::Type{T}, l::Integer, m::Integer) where T<:Real
    fac1 = one(T)
    for ii in (l-m+1):(l+m)
        fac1 *= convert(T, ii)
    end
    num = 2*convert(T, l) + 1
    den = 4*convert(T, π)
    return sqrt( num * inv(den) * inv(fac1) )
end
Nlm(l::Integer, m::Integer) = Nlm(Float64, l, m)

end # module Legendre
