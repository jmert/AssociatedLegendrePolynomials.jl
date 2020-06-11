# Implements a pre-computed table of recursion coefficients

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
            μ[m+1] = m == 0 ? zero(T) : coeff_μ(N(), T, m)
            # N.B.: Need access to mmax+1 but will never use m==0 term, so offset storage
            ν[m+1] = coeff_ν(N(), T, m+1)

            for l in (m+1):lmax
                α[l+1,m+1] = coeff_α(N(), T, l, m)
                β[l+1,m+1] = coeff_β(N(), T, l, m)
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

Base.convert(::Type{LegendreNormCoeff{N,T}}, norm::LegendreNormCoeff{N,T}) where {N,T} =
    norm
Base.convert(::Type{LegendreNormCoeff{N,T}}, norm::LegendreNormCoeff{N,S}) where {N,T,S} =
    LegendreNormCoeff{N,T}(norm)

Base.eltype(::Type{LegendreNormCoeff{N,T}}) where {N,T} = T

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

# Implements the legendre interface for the normalization

@inline function
initcond(::LegendreNormCoeff{N}, ::Type{T}) where {N<:AbstractLegendreNorm, T}
    return initcond(N(), T)
end

@propagate_inbounds function
coeff_μ(norm::LegendreNormCoeff, ::Type{T}, l::Integer) where T
    return norm.μ[l+1]
end

@propagate_inbounds function
coeff_ν(norm::LegendreNormCoeff, ::Type{T}, l::Integer) where T
    # N.B.: Storage is offset by 1 compared to other arrays. See constructor.
    return norm.ν[l]
end

@propagate_inbounds function
coeff_α(norm::LegendreNormCoeff, ::Type{T}, l::Integer, m::Integer) where T
    return norm.α[l+1,m+1]
end

@propagate_inbounds function
coeff_β(norm::LegendreNormCoeff, ::Type{T}, l::Integer, m::Integer) where T
    return norm.β[l+1,m+1]
end

function boundscheck_hook(norm::LegendreNormCoeff, lmax, mmax)
    @noinline _throw(α, lmax, mmax) = throw(BoundsError(α, (lmax+1, mmax+1)))
    lmax′,mmax′ = size(norm.α)
    (lmax < lmax′ && mmax < mmax′) || _throw(norm.α, lmax, mmax)
    nothing
end
