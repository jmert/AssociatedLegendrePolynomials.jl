import Base: checkindex, checkbounds_indices, DimOrInd, OneTo, Slice

@inline _similar(A::AbstractArray, ::Type{T}=eltype(A)) where {T} = similar(A, T, axes(A))
@inline _similar(A::Number, ::Type{T}=typeof(A)) where {T}        = Scalar{T}()

struct Work{T,N<:AbstractLegendreNorm,V<:AbstractArray{T}}
    norm::N
    z::V
    y¹::V
    y²::V
    pmm2::V
    pm::V
    plm2::V
    plm1::V
    pl::V
end
function Work(norm::AbstractLegendreNorm, Λ, x)
    T = promote_type(eltype(norm), eltype(Λ), eltype(x))
    z    = _similar(x, T)
    y¹   = _similar(x, T)
    y²   = _similar(x, T)
    pmm2 = _similar(x, T)
    pm   = _similar(x, T)
    plm2 = _similar(x, T)
    plm1 = _similar(x, T)
    pl   = _similar(x, T)
    return Work(norm, z, y¹, y², pmm2, pm, plm2, plm1, pl)
end

function _chkdomain(lmax, mmax)
    @noinline _chkdomain_throw_lmax(l) = throw(DomainError(l, "degree lmax must be non-negative"))
    @noinline _chkdomain_throw_mmax(m) = throw(DomainError(m, "order mmax must be non-negative and less than or equal to lmax"))

    0 ≤ lmax || _chkdomain_throw_lmax(lmax)
    0 ≤ mmax ≤ lmax || _chkdomain_throw_mmax(mmax)
    nothing
end

function _chkbounds(Λ, lmax, mmax, x)
    @noinline _chkbounds_throw_dims(M, N) = throw(DimensionMismatch(
        "Output has $M dimensions, expected $N to $(N+2)"))
    @noinline _chkbounds_throw_axes(Λ, x) = begin
        throw(DimensionMismatch(
            "Output has leading axes $(ntuple(i -> axes(Λ,i), ndims(x))), expected $(axes(x))"))
    end
    @noinline _chkbounds_throw_lmax() = throw(DimensionMismatch(
            "lmax incompatible with output array axes"))
    @noinline _chkbounds_throw_mmax() = throw(DimensionMismatch(
            "mmax incompatible with output array axes"))

    M = ndims(Λ)
    N = ndims(x)
    (0 ≤ M - N ≤ 2) || _chkbounds_throw_dims(M, N)
    # Leading dimensions of Λ are storage for the same dimensions as x
    axesΛ = axes(Λ)
    if N > 0
        axΛ = ntuple(i -> axesΛ[i], N)
        axx = axes(x)
        checkbounds_indices(Bool, axΛ, axx) || _chkbounds_throw_axes(Λ, x)
    end
    # Trailing dimensions of Λ are storage for range of ell and m, as needed
    dimΛ = ntuple(i -> axesΛ[i+N], M-N)
    if M - N > 0
        checkindex(Bool, dimΛ[1], OneTo(lmax+1)) || _chkbounds_throw_lmax()
    end
    if M - N > 1
        checkindex(Bool, dimΛ[2], OneTo(mmax+1)) || _chkbounds_throw_mmax()
    end
    nothing
end

function unsafe_legendre!(workornorm::Union{AbstractLegendreNorm,Work}, Λ, lmax, mmax, x)
    if ndims(x) > 1
        M = ndims(Λ)
        N = ndims(x)
        S = prod(size(x))
        x′ = reshape(x, S)
        Λ′ = reshape(Λ, S, ntuple(i->size(Λ,i+N), M-N)...)
    else
        x′ = x
        Λ′ = Λ
    end

    if workornorm isa AbstractLegendreNorm
        work = Work(workornorm, Λ′, x′)
        @inbounds _legendre_impl!(work, Λ′, lmax, mmax, x′)
    else
        @inbounds _legendre_impl!(workornorm, Λ′, lmax, mmax, x′)
    end
    return Λ
end

# Similar to `Base.fma` but defined for `Complex` arguments such that if all arguments are
# on the real axis, the result of `_fma(x, y, z) === _fma(real(x), real(y), real(z)) + 0im`.
#
# !!! note
#     The total operation for generic complex arguments is not fused but rather
#     proceeds via a series of FMA instructions (which include intermediate
#     rounding) --- the fusing is only strictly true for the all-real-axis case.
#
# Following definitions should be complete, but we comment out most of them since they
# are unused.
#_fma(z::Complex, w::Complex, x::Complex) =
#    Complex(_fma(real(z), real(w), -_fma(imag(z), imag(w), -real(x))),
#            _fma(real(z), imag(w),  _fma(imag(z), real(w),  imag(x))))
_fma(z::Complex, w::Complex, x::Real) =
    Complex(_fma(real(z), real(w), -_fma(imag(z), imag(w), -x)),
            _fma(real(z), imag(w), imag(z) * real(w)))
#_fma(x::Real, z::Complex, y::Number) = _fma(z, x, y)
#_fma(z::Complex, x::Real, y::Real) = Complex(_fma(real(z),x,y), imag(z)*x)
#_fma(z::Complex, x::Real, w::Complex) =
#    Complex(_fma(real(z),x,real(w)), _fma(imag(z),x,imag(w)))
#_fma(x::Real, y::Real, z::Complex) = Complex(_fma(x,y,real(z)), imag(z))
# Then cascade back to Base.fma
_fma(x, y, z) = Base.fma(x, y, z)

@propagate_inbounds function _legendre_impl!(work::Work, Λ, lmax, mmax, x)
    norm = work.norm
    z    = work.z
    y¹   = work.y¹
    y²   = work.y²
    pmm2 = work.pmm2
    pm   = work.pm
    plm2 = work.plm2
    plm1 = work.plm1
    pl   = work.pl

    local μ₁, μ₂
    T = eltype(z)
    M = ndims(x)
    N = ndims(Λ) - M
    I = CartesianIndices(map(Slice, axes(x)))

    @simd for ii in I
        z[ii]  = convert(T, x[ii])
        y²[ii] = -_fma(z[ii], z[ii], -one(real(T)))
        y¹[ii] = sqrt(y²[ii])
    end

    fill!(pm, initcond(norm, T))
    for m in 0:mmax
        @simd for ii in I
            if N == 2
                Λ[ii,m+1,m+1] = pm[ii]
            elseif N == 1 && m == mmax
                Λ[ii,m+1] = pm[ii]
            elseif N == 0 && m == lmax # == mmax
                Λ[ii] = pm[ii]
            end
        end
        m == lmax && break # exit for lmax == mmax == m

        if N == 2 || m == mmax
            # 1-term recurrence relation taking (l-1,l-1) -> (l,l-1) where l == m == m
            l = m + 1
            ν = coeff_ν(norm, real(T), l)
            @simd for ii in I
                plm1[ii] = pm[ii]
                pl[ii]   = ν * z[ii] * plm1[ii]
                if N == 2
                    Λ[ii,l+1,m+1] = pl[ii]
                elseif N == 1
                    Λ[ii,l+1] = pl[ii]
                elseif N == 0 && l == lmax
                    Λ[ii] = pl[ii]
                end
            end
            # no exit since following loop will be skipped for
            #   lmax == mmax + 1

            for l in m+2:lmax
                plm2, plm1, pl = plm1, pl, plm2
                # 2-term recurrence relation taking (l-1,m) -> (l, m)
                α = coeff_α(norm, real(T), l, m)
                β = coeff_β(norm, real(T), l, m)
                @simd for ii in I
                    pl[ii] = α * z[ii] * plm1[ii] - β * plm2[ii]
                    if N == 2
                        Λ[ii,l+1,m+1] = pl[ii]
                    elseif N == 1
                        Λ[ii,l+1] = pl[ii]
                    elseif N == 0 && l == lmax
                        Λ[ii] = pl[ii]
                    end
                end
            end
        end
        m == mmax && break # skip calculating (m,m) where m == mmax + 1

        # 1-term recurrence relation taking (m,m) -> (m+1,m+1)
        if iseven(m)
            pmm2, pm = pm, pmm2
            # Takes even m-2 to odd m-1 for following iteration where m will be odd
            μ₁ = coeff_μ(norm, real(T), m+1)
            @simd for ii in I
                pm[ii] = -μ₁ * y¹[ii] * pmm2[ii]
            end
        else
            # Takes even m-2 to even m for following iteration where m will be even again
            μ₂ = coeff_μ(norm, real(T), m+1)
            @simd for ii in I
                pm[ii] = μ₁ * μ₂ * y²[ii] * pmm2[ii]
            end
        end
    end

    return Λ
end

"""
    legendre!(norm::AbstractLegendreNorm, Λ, l::Integer, m::Integer, x)

Fills the array `Λ` with the Legendre polynomial values ``N_ℓ^m P_ℓ^m(x)``
up to/of degree(s) `l` and order(s) `m` for the normalization scheme `norm`.
`Λ` must be an array with between 0 and 2 more dimensions than `x`, with the leading
dimensions having the same shape as `x`.

- If `ndims(Λ) == ndims(x)`, then `Λ` is filled with the polynomial values at `x` for
  degree `l` and order `m`.
- If `ndims(Λ) == ndims(x) + 1`, then `l` is interpreted as `lmax`, and `Λ` filled with
  polynomial values for all degrees `0 ≤ l ≤ lmax` of order `m`.
- If `ndims(Λ) == ndims(x) + 2`, then `l` is interpreted as `lmax` and `m` as `mmax`,
  and `Λ` is filled with polynomial values for all degrees `0 ≤ l ≤ lmax` and orders
  `0 ≤ m ≤ min(mmax, l)`.
"""
function legendre!(norm::AbstractLegendreNorm, Λ, l::Integer, m::Integer, x)
    _chkdomain(l, m)
    boundscheck_hook(norm, l, m)
    _chkbounds(Λ, l, m, x)
    return unsafe_legendre!(norm, Λ, l, m, x)
end

"""
    p = legendre(norm::AbstractLegendreNorm, l::DimOrInd, m::DimOrInd, x)

Computes the associated Legendre polynomials ``N_ℓ^m P_ℓ^m(x)`` of degree(s) `l` and
order(s) `m` at `x` for the normalization scheme `norm`.

- If `l` and `m` are integers, returns a single value.
- If `l` is a range `0:lmax` and `m` an integer, returns the vector of values for order
  `m` and all degrees `0 ≤ l ≤ lmax`.
- If `l` is a range `0:lmax` and `m` is a range `0:mmax`, returns the matrix of values
  for all degrees `0 ≤ l ≤ lmax` and orders `0 ≤ m ≤ mmax`.

Note that in both the second and third cases, the ranges must have a first index of 0.
"""
function legendre end

# Scalar argument, scalar output case handled separately.
function legendre(norm::AbstractLegendreNorm, l::Integer, m::Integer, x::Number)
    _chkdomain(l, m)
    boundscheck_hook(norm, l, m)
    Λ = _similar(x)
    unsafe_legendre!(norm, Λ, l, m, x)
    return Λ[]
end

function legendre(norm::AbstractLegendreNorm, l::DimOrInd, m::DimOrInd, x)
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

    Λ = zeros(eltype(x), axes(x)..., size(l)..., size(m)...)
    unsafe_legendre!(norm, Λ, lmax, mmax, x)
    return Λ
end

# Make normalizations callable with similar syntax as the legendre[!] functions
@static if VERSION < v"1.3.0-alpha"
    # Prior to julia-1.3.0, abstract types could not have methods attached to them.
    # Provide for just the defined normalization types.
    for N in (LegendreUnitNorm,LegendreSphereNorm,LegendreOrthoNorm,LegendreFourPiNorm,LegendreNormCoeff)
        @eval begin
            @inline function (norm::$N)(l, m, x)
                return legendre(norm, l, m, x)
            end
            @inline function (norm::$N)(Λ, l::Integer, m::Integer, x)
                return legendre!(norm, Λ, l, m, x)
            end
        end
    end
else
    @inline function (norm::AbstractLegendreNorm)(l, m, x)
        return legendre(norm, l, m, x)
    end
    @inline function (norm::AbstractLegendreNorm)(Λ, l::Integer, m::Integer, x)
        return legendre!(norm, Λ, l, m, x)
    end
end
