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

import Base:
    @boundscheck, @propagate_inbounds,
    broadcastable, convert, eltype

module Bcast
    import Base: isassigned, getindex, setindex!, copyto!, @propagate_inbounds
    import Base.Broadcast: Broadcasted, AbstractArrayStyle, dotview, materialize

    """
        BScalar{T} <: Ref{T}

    A scalar storage location within broadcast expressions, very similar to `Base.RefValue`,
    but where assignment in a broadcasted context (for dimensionally-compatible expressions)
    is supported. For example:
    ```
    f!(b::BScalar, x) = (b .= sin.(x); b)
    b = BScalar{Float64}();
    f!(b, π/2)
    ```
    """
    struct BScalar{T} <: Ref{T}
        x::Base.RefValue{T}
        BScalar{T}() where {T} = new(Ref{T}())
        BScalar{T}(x) where {T} = new(x)
    end
    BScalar(x::T) where {T} = BScalar{T}(x)
    # Minimal copy of RefValue's methods
    @propagate_inbounds isassigned(x::BScalar) = isassigned(x.x)
    @propagate_inbounds getindex(b::BScalar) = b.x[]
    @propagate_inbounds setindex!(b::BScalar, x) = (b.x[] = x; b)
    # Broadcasting extensions
    @propagate_inbounds dotview(A::BScalar, ::CartesianIndex{0}) = A
    @propagate_inbounds dotview(A::BScalar, ::CartesianIndices{0,Tuple{}}) = A
    @propagate_inbounds copyto!(dest::BScalar, bc::Broadcasted{<:AbstractArrayStyle{0}}) =
        dest[] = materialize(bc)
end # module Bcast
import .Bcast: BScalar

"""
    abstract type AbstractLegendreNorm end

Abstract supertype for normalization conditions of the Associated Legendre polynomials.

# Example
```jldoctest
julia> using InteractiveUtils; subtypes(AbstractLegendreNorm)
3-element Array{Any,1}:
 LegendreNormCoeff
 LegendreSphereNorm
 LegendreUnitNorm
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
LegendreNormCoeff{LegendreSphereNorm,Float64} for lmax = 1, mmax = 1 with coefficients:
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
        (lmax ≥ 0) || throw(DomainError(lmax, "lmax must be positive"))
        (0 ≤ mmax ≤ lmax) || throw(DomainError(mmax, "mmax must be bounded in 0 to $lmax"))

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

    function LegendreNormCoeff{N,T1}(norm::LegendreNormCoeff{N,T2}) where {N,T1,T2}
        return new(convert(Vector{T1}, norm.μ),
                   convert(Vector{T1}, norm.ν),
                   convert(Matrix{T1}, norm.α),
                   convert(Matrix{T1}, norm.β))
    end
end

LegendreNormCoeff{N,T}(lmax::Integer) where {N,T} = LegendreNormCoeff{N,T}(lmax, lmax)

convert(::Type{LegendreNormCoeff{N,T}}, norm::LegendreNormCoeff{N,T}) where {N,T} = norm
convert(::Type{LegendreNormCoeff{N,T}}, norm::LegendreNormCoeff{N,S}) where {N,T,S} =
        LegendreNormCoeff{N,T}(norm)

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

# Define element types as a way to conveniently promote, even for the trait types.
# Union{} will always lose out in type promotion.
eltype(::Type{<:AbstractLegendreNorm}) = Union{}
eltype(::Type{LegendreNormCoeff{N,T}}) where {N,T} = T

# Normalizations are scalars in a broadcasting context
broadcastable(x::AbstractLegendreNorm) = Ref(x)

# Improve printing somewhat
Base.show(io::IO, norm::LegendreNormCoeff{N,T}) where {N,T} =
    print(io, LegendreNormCoeff, "{$N,$T}")
function Base.show(io::IO, ::MIME"text/plain", N::LegendreNormCoeff)
    lmax,mmax = size(N.α) .- 1
    println(io, N, " for lmax = $lmax, mmax = $mmax with coefficients:")
    io′ = IOContext(io, :compact => true)
    println(io′, "    μ: ", N.μ)
    println(io′, "    ν: ", N.ν)
    println(io′, "    α: ", N.α)
    println(io′, "    β: ", N.β)
end

"""
    Plm_00(::N, ::Type{T}) where {N<:AbstractLegendreNorm, T}

Returns the initial condition ``P_0^0(x)`` for the associated Legendre recursions based
on the normalization choice `N` for numeric type `T`.
"""
function Plm_00 end

"""
    Plm_μ(norm::N, ::Type{T}, l::Integer) where {N<:AbstractLegendreNorm, T}

Returns the coefficient ``μ_ℓ`` for the single-term recursion relation
```math
    P_{ℓ+1}^{ℓ+1}(x) = -μ_{ℓ+1} \\sqrt{1-x^2} P_ℓ^ℓ(x)
```
where ``μ_ℓ`` is appropriate for the choice of normalization `N`.
"""
function Plm_μ end

"""
    Plm_ν(norm::N, ::Type{T}, l::Integer) where {N<:AbstractLegendreNorm, T}

Returns the coefficient ``ν_ℓ`` for the single-term recursion relation
```math
    P_{ℓ+1}^ℓ(x) = ν_{ℓ+1} x P_ℓ^ℓ(x)
```
where ``ν_ℓ`` is appropriate for the choice of normalization `N`.
"""
function Plm_ν end

"""
    Plm_α(norm::N, ::Type{T}, l::Integer, m::Integer) where {N<:AbstractLegendreNorm, T}

Returns the coefficient ``α_ℓ^m`` for the two-term recursion relation
```math
    P_{ℓ+1}^{m}(x) = α_{ℓ+1}^m x P_ℓ^m(x) - β_{ℓ+1}^m P_{ℓ-1}^m(x)
```
where ``α_ℓ^m`` is appropriate for the choice of normalization `N`.
"""
function Plm_α end

"""
    Plm_β(norm::N, ::Type{T}, l::Integer, m::Integer) where {N<:AbstractLegendreNorm, T}

Returns the coefficient ``β_ℓ^m`` for the two-term recursion relation
```math
    P_{ℓ+1}^{m}(x) = α_{ℓ+1}^m x P_ℓ^m(x) - β_{ℓ+1}^m P_{ℓ-1}^m(x)
```
where ``β_ℓ^m`` is appropriate for the choice of normalization `N`.
"""
function Plm_β end

@inline function
Plm_00(::LegendreUnitNorm, ::Type{T}) where T
    return one(T)
end

@inline function
Plm_00(::LegendreSphereNorm, ::Type{T}) where T
    # comparing this against
    #   convert(T, inv(sqrt(4*convert(BigFloat, π))))
    # shows that this is exact within Float64 precision
    return inv(sqrt(4 * convert(T, π)))
end

@inline function
Plm_00(::LegendreNormCoeff{N}, ::Type{T}) where {N<:AbstractLegendreNorm, T}
    return Plm_00(N(), T)
end

@inline function
Plm_μ(::LegendreUnitNorm, ::Type{T}, m::Integer) where T
    return convert(T, 2m - 1)
end

@inline function
Plm_μ(::LegendreSphereNorm, ::Type{T}, m::Integer) where T
    return sqrt(one(T) + inv(convert(T, 2m)))
end

@propagate_inbounds function
Plm_μ(norm::LegendreNormCoeff, ::Type{T}, m::Integer) where T
    return norm.μ[m+1]
end

@inline function
Plm_ν(::LegendreUnitNorm, ::Type{T}, m::Integer) where T
    return convert(T, 2m + 1)
end

@inline function
Plm_ν(::LegendreSphereNorm, ::Type{T}, m::Integer) where T
    return sqrt(convert(T, 2m + 3))
end

@propagate_inbounds function
Plm_ν(norm::LegendreNormCoeff, ::Type{T}, m::Integer) where T
    return norm.ν[m+1]
end

@inline function
Plm_α(::LegendreUnitNorm, ::Type{T}, l::Integer, m::Integer) where T
    return convert(T, 2l-1) * inv(convert(T, l-m))
end

@inline function
Plm_α(::LegendreSphereNorm, ::Type{T}, l::Integer, m::Integer) where T
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
Plm_α(norm::LegendreNormCoeff, ::Type{T}, l::Integer, m::Integer) where T
    return norm.α[l+1,m+1]
end

@inline function
Plm_β(::LegendreUnitNorm, ::Type{T}, l::Integer, m::Integer) where T
    return convert(T, l+m-1) * inv(convert(T, l-m))
end

@inline function
Plm_β(::LegendreSphereNorm, ::Type{T}, l::Integer, m::Integer) where T
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
Plm_β(norm::LegendreNormCoeff, ::Type{T}, l::Integer, m::Integer) where T
    return norm.β[l+1,m+1]
end

# Named recursion relations

@propagate_inbounds function
_1term_raise_lm(norm::AbstractLegendreNorm, m::Integer, x::T, y::T, plm::T) where T
    μ = Plm_μ(norm, T, m+1)
    return -μ * y * plm
end
@propagate_inbounds function
_1term_raise_l(norm::AbstractLegendreNorm, m::Integer, x::T, plm::T) where T
    ν = Plm_ν(norm, T, m)
    return ν * x * plm
end
@propagate_inbounds function
_2term_raise_l(norm::AbstractLegendreNorm, l::Integer, m::Integer, x::T,
               plm::T, plm1m::T) where T
    α = Plm_α(norm, T, l+1, m)
    β = Plm_β(norm, T, l+1, m)
    return α * x * plm - β * plm1m
end

@noinline function _chkdomain(lmax, mmax)
    0 ≤ lmax || throw(DomainError(lmax, "degree lmax must be non-negative"))
    0 ≤ mmax ≤ lmax || throw(DomainError(mmax,
            "order mmax must be non-negative and less than lmax"))
    nothing
end
@noinline function _chkbounds(Λ, lmax, mmax, x)
    M = ndims(Λ)
    N = ndims(x)
    # Leading dimensions of Λ are storage for the same dimensions as x
    if N > 0
        szΛ = ntuple(i -> size(Λ, i), N)
        szx = size(x)
        all(szΛ .>= szx) || throw(DimensionMismatch(
                "Output storage has leading dimensions of size $szΛ, need at least $szx"))
    end
    # Trailing dimensions of Λ are storage for range of ell and m, as needed
    dimΛ = ntuple(i -> size(Λ, i+N), M-N)
    if length(dimΛ) > 0
        lmax < dimΛ[1] || throw(DimensionMismatch(
                "lmax incompatible with output array size"))
    end
    if length(dimΛ) > 1
        mmax < dimΛ[2] || throw(DimensionMismatch(
                "mmax incompatible with output array size"))
    end
    nothing
end

@propagate_inbounds function _legendre!(norm, Λ, lmax, mmax, x)
    @boundscheck _chkdomain(lmax, mmax)
    @boundscheck _chkbounds(Λ, lmax, mmax, x)
    @inbounds _legendre_impl!(norm, Λ, lmax, mmax, x)
end

@inline _similar(A::AbstractArray) = similar(A, size(A))
@inline _similar(A::BScalar)       = BScalar{eltype(A)}()
@inline _similar(A::Number)        = BScalar{typeof(A)}()
@propagate_inbounds function _legendre_impl!(norm::AbstractLegendreNorm, Λ, lmax, mmax, x)
    TΛ = eltype(Λ)
    TV = eltype(x)
    T = promote_type(eltype(norm), TΛ, TV)

    z    = _similar(x)
    y    = _similar(x)
    pmp1 = _similar(x)
    pm   = _similar(x)
    plm1 = _similar(x)
    pl   = _similar(x)
    plp1 = _similar(x)

    z .= convert.(T, x)
    y .= @fastmath sqrt.(one.(T) .- z .^ 2)

    M = ndims(x)
    N = ndims(Λ) - M
    sz = size(x)
    I = CartesianIndices(sz)

    pmp1 .= Plm_00(norm, T)
    for m in 0:mmax
        pm, pmp1 = pmp1, pm
        if m < mmax
            pmp1 .= _1term_raise_lm.(norm, m, z, y, pm)
        end

        if N == 2
            Λ[I,m+1,m+1] .= pm
        elseif N == 1
            Λ[I,m+1] .= (m != mmax) ? zero(T) : pm
            m != mmax && continue
        elseif N == 0
            m != mmax && continue
            if lmax == m
                Λ[I] .= pm
                return Λ
            end
        end

        pl   .= pm
        plp1 .= _1term_raise_l.(norm, m, z, pm)
        for l in m+1:lmax
            plm1, pl, plp1 = pl, plp1, plm1
            if l < lmax
                plp1 .= _2term_raise_l.(norm, l, m, z, pl, plm1)
            end

            if N == 2
                Λ[I,l+1,m+1] .= pl
            elseif N == 1
                Λ[I,l+1] .= pl
            elseif N == 0 && m == mmax && l == lmax
                Λ[I] .= pl
            end
        end
    end

    return Λ
end

"""
    p = legendre(norm::AbstractLegendreNorm, l::Integer, x::Number)

Computes the scalar value ``p = N_ℓ P_ℓ(x)``, where ``P_ℓ(x)`` is the Legendre
polynomial of degree `l` at `x` and ``N_ℓ`` is the normalization scheme `norm`.
"""
@propagate_inbounds function legendre(
        norm::AbstractLegendreNorm, l::Integer, x::Number)
    Λ = BScalar{typeof(x)}()
    _legendre!(norm, Λ, l, 0, x)
    return Λ[]
end

"""
    p = legendre(norm::AbstractLegendreNorm, l::Integer, m::Integer, x::Number)

Computes the scalar value ``p = N_ℓ^m P_ℓ^m(x)``, where ``P_ℓ^m(x)`` is the associated
Legendre polynomial of degree `l` and order `m` at `x` and ``N_ℓ^m`` the normalization
scheme `norm`.
"""
@propagate_inbounds function legendre(
        norm::AbstractLegendreNorm, l::Integer, m::Integer, x::Number)
    Λ = BScalar{typeof(x)}()
    _legendre!(norm, Λ, l, m, x)
    return Λ[]
end

"""
    legendre!(norm::AbstractLegendreNorm, Λ, l::Integer, m::Integer, x)

Fills the array `Λ` with the Legendre polynomial values ``N_ℓ^m P_ℓ^m(x)``, where ``N_ℓ^m``
is the normalization scheme `norm`.
"""
@propagate_inbounds function legendre!(
        norm::AbstractLegendreNorm,
        Λ, l::Integer, m::Integer, x)
    return _legendre!(norm, Λ, l, m, x)
end

# Make coefficient cache objects callable with similar syntax as the legendre[!] functions
@propagate_inbounds function (norm::LegendreNormCoeff)(l::Integer, x)
    return legendre(norm, l, x)
end
@propagate_inbounds function (norm::LegendreNormCoeff)(l::Integer, m::Integer, x)
    return legendre(norm, l, m, x)
end
@propagate_inbounds function (norm::LegendreNormCoeff)(Λ, m::Integer, x)
    lmax = size(norm.α,1) - 1
    return legendre!(norm, Λ, lmax, m, x)
end
@propagate_inbounds function (norm::LegendreNormCoeff)(Λ, x)
    lmax,mmax = size(norm.α) .- 1
    return legendre!(norm, Λ, lmax, mmax, x)
end


"""
    p = Pl(l::Integer, x)

Computes the scalar value ``p = P_ℓ(x)``, where ``P_ℓ(x)`` is the Legendre polynomial of
degree `l` at `x`.
"""
@inline Pl(l::Integer, x) = legendre(LegendreUnitNorm(), l, x)

"""
    p = Plm(l::Integer, m::Integer, x)

Computes the scalar value ``p = P_ℓ^m(x)``, where ``P_ℓ^m(x)`` is the associated Legendre
polynomial of degree `l` and order `m` at `x`.
"""
@inline Plm(l::Integer, m::Integer, x) = legendre(LegendreUnitNorm(), l, m, x)

"""
    λ = λlm(l::Integer, m::Integer, x)

Computes the scalar value ``λ = λ_ℓ^m(x)``, where ``λ_ℓ^m(x)`` is the spherical-harmonic
normalized associated Legendre polynomial of degree `l` and order `m` at `x`.
"""
@inline λlm(l::Integer, m::Integer, x) = legendre(LegendreSphereNorm(), l, m, x)

"""
    Pl!(P, l::Integer, x)

Fills the vector `P` with the Legendre polynomial values ``P_ℓ(x)`` for all degrees
`0 ≤ ℓ ≤ lmax` at `x`.
"""
@inline Pl!(P, l::Integer, x) =
    legendre!(LegendreUnitNorm(), P, l, 0, x)

"""
    Plm!(P, l::Integer, m::Integer, x)

Fills the vector `P` with the Legendre polynomial values ``P_ℓ^m(x)`` for all degrees
`0 ≤ ℓ ≤ lmax` and constant order `m` at `x`.
"""
@inline Plm!(P, l::Integer, m::Integer, x) =
    legendre!(LegendreUnitNorm(), P, l, m, x)


"""
    λlm!(Λ, l::Integer, m::Integer, x)

Fills the vector `Λ` with the spherical harmonic normalized associated Legendre polynomial
values ``λ_ℓ^m(x)`` for all degrees `0 ≤ ℓ ≤ lmax` and constant order `m` at `x`.
"""
@inline λlm!(Λ, l::Integer, m::Integer, x) =
    legendre!(LegendreSphereNorm(), Λ, l, m, x)

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

end # module Legendre
