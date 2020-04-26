@inline _chkdomainnorm(norm::AbstractLegendreNorm, lmax, mmax) = nothing
@noinline function _chkdomainnorm(norm::LegendreNormCoeff, lmax, mmax)
    lmax′,mmax′ = size(norm.α)
    (lmax < lmax′ && mmax < mmax′) || throw(BoundsError(norm.α, (lmax+1, mmax+1)))
end
@noinline function _chkdomain(lmax, mmax)
    0 ≤ lmax || throw(DomainError(lmax, "degree lmax must be non-negative"))
    0 ≤ mmax ≤ lmax || throw(DomainError(mmax,
            "order mmax must be non-negative and less than or equal to lmax"))
    nothing
end
@noinline function _chkbounds(Λ, lmax, mmax, x)
    M = ndims(Λ)
    N = ndims(x)
    if !(0 ≤ M - N ≤ 2)
        throw(DimensionMismatch(
                "Output storage has $M dimensions; expected $N to $(N+2)"))
    end
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
    @boundscheck _chkdomainnorm(norm, lmax, mmax)
    @boundscheck _chkbounds(Λ, lmax, mmax, x)
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

@inline _similar(A::AbstractArray) = similar(A, size(A))
@inline _similar(A::Number)        = Scalar{typeof(A)}()
@propagate_inbounds function _legendre_impl!(norm::AbstractLegendreNorm, Λ, lmax, mmax, x)
    TΛ = eltype(Λ)
    TV = eltype(x)
    T = promote_type(eltype(norm), TΛ, TV)

    z    = _similar(x)
    y¹   = _similar(x)
    y²   = _similar(x)
    pm   = _similar(x)
    pmp2 = _similar(x)
    plm1 = _similar(x)
    pl   = _similar(x)
    plp1 = _similar(x)

    @. z = convert(T, x)
    @. y² = -fma(z, z, -one(T))
    @. y¹ = sqrt(y²)

    M = ndims(x)
    N = ndims(Λ) - M
    sz = size(x)
    I = CartesianIndices(sz)
    mmax′ = mmax - mod(unsigned(mmax), 2)

    local μ₁, μ₂
    pmp2 .= Plm_00(norm, T)
    for m in 0:2:mmax′
        pm, pmp2 = pmp2, pm
        pmp1 = pmp2 # name alias

        # First iteration: takes even m to m+1
        if N == 2
            Λ[I,m+1,m+1] = pm
        elseif N == 1 && m == mmax
            Λ[I,m+1] = pm
        elseif N == 0 && lmax == m
            Λ[I] = pm
            return Λ
        end
        if N == 2 || m == mmax
            # 1-term recurrence relation taking (m,m) -> (m,m+1)
            ν = Plm_ν(norm, T, m)
            @. pl   = pm
            @. plp1 = ν * z * pl
            for l in m+1:lmax
                plm1, pl, plp1 = pl, plp1, plm1
                if l < lmax
                    # 2-term recurrence relation taking (l,m) -> (l+1, m)
                    α = Plm_α(norm, T, l+1, m)
                    β = Plm_β(norm, T, l+1, m)
                    @. plp1 = α * z * pl - β * plm1
                end

                if N == 2
                    Λ[I,l+1,m+1] = pl
                elseif N == 1
                    Λ[I,l+1] = pl
                elseif N == 0 && l == lmax
                    Λ[I] = pl
                    return
                end
            end
        end
        if m < mmax
            # 1-term recurrence relation taking (m,m) -> (m+1,m+1)
            μ₁ = Plm_μ(norm, T, m+1)
            @. pmp1 = -μ₁ * y¹ * pm
        end

        m == mmax == mmax′ && break

        # Second iteration: takes even m to m+2

        if N == 2
            Λ[I,m+2,m+2] = pmp1
        elseif N == 1 && m+1 == mmax
            Λ[I,m+2] = pmp1
        elseif N == 0 && m+1 == lmax
            Λ[I] = pmp1
            return Λ
        end
        if N == 2 || m+1 == mmax
            # 1-term recurrence relation taking (m+1,m+1) -> (m+1,m+2)
            ν = Plm_ν(norm, T, m+1)
            @. pl   = pmp1
            @. plp1 = ν * z * pl
            for l in m+2:lmax
                plm1, pl, plp1 = pl, plp1, plm1
                if l < lmax
                    # 2-term recurrence relation taking (l,m+1) -> (l+1, m+1)
                    α = Plm_α(norm, T, l+1, m+1)
                    β = Plm_β(norm, T, l+1, m+1)
                    @. plp1 = α * z * pl - β * plm1
                end

                if N == 2
                    Λ[I,l+1,m+2] = pl
                elseif N == 1
                    Λ[I,l+1] = pl
                elseif N == 0 && l == lmax
                    Λ[I] = pl
                    return
                end
            end
        end
        if m+1 < mmax
            # 1-term recurrence relation applied twice taking (m,m) -> (m+2,m+2)
            μ₂ = Plm_μ(norm, T, m+2)
            @. pmp2 = μ₁ * μ₂ * y² * pm
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
    _chkdomainnorm(norm, l, m)
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
@propagate_inbounds function legendre!(
        norm::AbstractLegendreNorm,
        Λ, l::Integer, m::Integer, x)
    return _legendre!(norm, Λ, l, m, x)
end
