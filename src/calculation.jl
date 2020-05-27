import Base: checkindex, checkbounds_indices, OneTo, Slice

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

@propagate_inbounds function _legendre!(norm, Λ, lmax, mmax, x)
    @boundscheck _chkdomain(lmax, mmax)
    @boundscheck _chkbounds(Λ, lmax, mmax, x)
    @boundscheck boundscheck_hook(norm, lmax, mmax)
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
    @inbounds _legendre_impl!(norm, Λ′, lmax, mmax, x′)
    return Λ
end

@inline _similar(A::AbstractArray, ::Type{T}=eltype(A)) where {T} = similar(A, T, axes(A))
@inline _similar(A::Number, ::Type{T}=typeof(A)) where {T}        = Scalar{T}()
@propagate_inbounds function _legendre_impl!(norm::AbstractLegendreNorm, Λ, lmax, mmax, x)
    TΛ = eltype(Λ)
    TV = eltype(x)
    T = promote_type(eltype(norm), TΛ, TV)

    z    = _similar(x, T)
    y¹   = _similar(x, T)
    y²   = _similar(x, T)
    pmm2 = _similar(x, T)
    pm   = _similar(x, T)
    plm2 = _similar(x, T)
    plm1 = _similar(x, T)
    pl   = _similar(x, T)

    @. z = convert(T, x)
    @. y² = -fma(z, z, -one(T))
    @. y¹ = sqrt(y²)

    M = ndims(x)
    N = ndims(Λ) - M
    I = CartesianIndices(map(Slice, axes(x)))

    local μ₁, μ₂
    pm .= Plm_00(norm, T)
    for m in 0:mmax
        if N == 2
            Λ[I,m+1,m+1] = pm
        elseif N == 1 && m == mmax
            Λ[I,m+1] = pm
        elseif N == 0 && m == lmax # == mmax
            Λ[I] = pm
        end
        m == lmax && break # exit for lmax == mmax == m

        if N == 2 || m == mmax
            # 1-term recurrence relation taking (l-1,l-1) -> (l,l-1) where l == m == m
            l = m + 1
            ν = Plm_ν(norm, T, l)
            @. plm1 = pm
            @. pl   = ν * z * plm1
            if N == 2
                Λ[I,l+1,m+1] = pl
            elseif N == 1
                Λ[I,l+1] = pl
            elseif N == 0 && l == lmax
                Λ[I] = pl
            end
            # no exit since following loop will be skipped for
            #   lmax == mmax + 1

            for l in m+2:lmax
                plm2, plm1, pl = plm1, pl, plm2
                # 2-term recurrence relation taking (l-1,m) -> (l, m)
                α = Plm_α(norm, T, l, m)
                β = Plm_β(norm, T, l, m)
                @. pl = α * z * plm1 - β * plm2

                if N == 2
                    Λ[I,l+1,m+1] = pl
                elseif N == 1
                    Λ[I,l+1] = pl
                elseif N == 0 && l == lmax
                    Λ[I] = pl
                end
            end
        end
        m == mmax && break # skip calculating (m,m) where m == mmax + 1

        # 1-term recurrence relation taking (m,m) -> (m+1,m+1)
        if iseven(m)
            pmm2, pm = pm, pmm2
            # Takes even m-2 to odd m-1 for following iteration where m will be odd
            μ₁ = Plm_μ(norm, T, m+1)
            @. pm = -μ₁ * y¹ * pmm2
        else
            # Takes even m-2 to even m for following iteration where m will be even again
            μ₂ = Plm_μ(norm, T, m+1)
            @. pm = μ₁ * μ₂ * y² * pmm2
        end
    end

    return Λ
end

"""
    p = legendre(norm::AbstractLegendreNorm, l::Integer, m::Integer, x::Number)
    P = legendre.(norm::AbstractLegendreNorm, l, m, x)

Computes the associated Legendre polynomial ``N_ℓ^m P_ℓ^m(x)`` of degree `l` and
order `m` at `x` for the normalization scheme `norm`.

With broadcasting syntax, the polynomials can be computed over any iterable `x`.
Furthermore,
- If `l isa Integer && m isa Integer`, then the output `P` has the same shape as `x` and
  is filled with the polynomial values of order `l` and degree `m`.
- If `l isa UnitRange && m isa Integer`, then `l` is interpreted as `lmax`, and the output
  `P` has one more dimension than `x` with the trailing dimension spanning the degrees
  `0 ≤ l ≤ lmax`.
- If `l isa UnitRange && m isa UnitRange`, then `l` is interpreted as `lmax` and `m` as
  `mmax`, and the output `P` has two more dimensions than `x` with the trailing dimensions
  spanning the degrees `0 ≤ l ≤ lmax` and orders `0 ≤ m ≤ mmax`, respectively.

Note that in second and third case, the `UnitRange`s must satisify `first(l) == 0` and
`first(m) == 0`.
"""
function legendre(norm::AbstractLegendreNorm, l::Integer, m::Integer, x::Number)
    Λ = _similar(x)
    _chkdomain(l, m)
    boundscheck_hook(norm, l, m)
    @inbounds _legendre!(norm, Λ, l, m, x)
    return Λ[]
end

"""
    legendre!(norm::AbstractLegendreNorm, Λ, l::Integer, m::Integer, x)

Fills the array `Λ` with the Legendre polynomial values ``N_ℓ^m P_ℓ^m(x)``
up to/of degree `l` and order `m` for the normalization scheme `norm`.
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
function legendre!(
        norm::AbstractLegendreNorm,
        Λ, l::Integer, m::Integer, x)
    return _legendre!(norm, Λ, l, m, x)
end

# Make normalizations callable with similar syntax as the legendre[!] functions
@static if VERSION < v"1.3.0-alpha"
    # Prior to julia-1.3.0, abstract types could not have methods attached to them.
    # Provide for just the defined normalization types.
    for N in (LegendreUnitNorm,LegendreSphereNorm,LegendreNormCoeff)
        @eval begin
            @inline function (norm::$N)(l::Integer, m::Integer, x)
                return legendre(norm, l, m, x)
            end
            @inline function (norm::$N)(Λ, l::Integer, m::Integer, x)
                return legendre!(norm, Λ, l, m, x)
            end
        end
    end
else
    @inline function (norm::AbstractLegendreNorm)(l::Integer, m::Integer, x)
        return legendre(norm, l, m, x)
    end
    @inline function (norm::AbstractLegendreNorm)(Λ, l::Integer, m::Integer, x)
        return legendre!(norm, Λ, l, m, x)
    end
end
