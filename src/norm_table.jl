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
Plm_00(::LegendreNormCoeff{N}, ::Type{T}) where {N<:AbstractLegendreNorm, T}
    return Plm_00(N(), T)
end

@propagate_inbounds function
Plm_μ(norm::LegendreNormCoeff, ::Type{T}, m::Integer) where T
    return norm.μ[m+1]
end

@propagate_inbounds function
Plm_ν(norm::LegendreNormCoeff, ::Type{T}, m::Integer) where T
    return norm.ν[m+1]
end

@propagate_inbounds function
Plm_α(norm::LegendreNormCoeff, ::Type{T}, l::Integer, m::Integer) where T
    return norm.α[l+1,m+1]
end

@propagate_inbounds function
Plm_β(norm::LegendreNormCoeff, ::Type{T}, l::Integer, m::Integer) where T
    return norm.β[l+1,m+1]
end
