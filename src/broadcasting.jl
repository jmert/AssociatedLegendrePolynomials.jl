# Specialize broadcasting of all of the non-modifying interfaces
import Base.broadcasted

@inline broadcasted(norm::T, l, m, x) where {T<:AbstractLegendreNorm} =
    broadcasted(legendre, norm, l, m, x)

# N.B. The broadcasting functions are near replicas of legendre() functions in
#      calculation.jl. Done so that bounds checking can occur before any array allocations.

function broadcasted(::typeof(legendre),
        norm::AbstractLegendreNorm, l::Integer, m::Integer, x)
    _chkdomain(l, m)
    boundscheck_hook(norm, l, m)
    z = Broadcast.materialize(x)
    Λ = _similar(z)
    _legendre!(norm, Λ, l, m, z)
    return (ndims(Λ) == 0) ? Λ[] : Λ
end

function broadcasted(::typeof(legendre), norm::AbstractLegendreNorm, l::DimOrInd, m::DimOrInd, x)
    if l isa AbstractUnitRange
        first(l) == 0 || throw(ArgumentError("Range of degrees l must start at 0"))
    end
    if m isa AbstractUnitRange
        if !(l isa AbstractUnitRange)
            throw(ArgumentError("Range of orders m requires range of degrees l"))
        end
        first(m) == 0 || throw(ArgumentError("Range of orders m must start at 0"))
    end
    lmax = l isa AbstractUnitRange ? last(l) : l
    mmax = m isa AbstractUnitRange ? last(m) : m
    _chkdomain(lmax, mmax)
    boundscheck_hook(norm, lmax, mmax)

    z = Broadcast.materialize(x)
    Λ = zeros(eltype(z), axes(z)..., size(l)..., size(m)...)
    _legendre!(norm, Λ, lmax, mmax, z)
    return Λ
end
