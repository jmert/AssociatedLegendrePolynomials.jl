"""
Collections of functions which compute the associated Legendre functions.

Based on implementation described in Limpanuparb and Milthorpe (2014)
*“Associated Legendre Polynomials and Spherical Harmonics Computation for
Chemistry Applications”* arXiv:1410.1748v1
"""
module Legendre
    import Base.@boundscheck, Base.@propagate_inbounds

    export
        PlmPlan, PlmSphericalPlan,
        Pl, Pl!,
        Plms, Plm, Plm!

    #########################################################################
    # LEGENDRE POLYNOMIALS
    #########################################################################

    doc"""
        Pl(lmax, x)

    Returns the Legendre polynomials up to degree `lmax` evaluated at `x`.
    """
    function Pl end

    function Pl(lmax::Integer, x::Number)
        pl = Array{eltype(x)}(1, lmax+1)
        @inbounds Pl!(pl, lmax, x)
        return reshape(pl, lmax+1)
    end
    function Pl(lmax::Integer, x::Vector{T} where T<:Number)
        pl = Array{eltype(x)}(length(x), lmax+1)
        @inbounds Pl!(pl, lmax, x)
        return pl
    end

    @inline function _chkbounds_Pl(pl, lmax, x)
        lmax >= 0 || throw(DimensionMismatch("lmax must be non-negative"))
        size(pl) >= (length(x),lmax+1) || throw(DimensionMismatch())
    end

    doc"""
        Pl!(pl, lmax, x)

    Fills `pl` with the Legendre polynomials up to degree `lmax` evaluated at
    `x`.
    """
    @propagate_inbounds function Pl!(pl::Matrix{T1} where T1<:Number, lmax::Integer, x::Union{T2,Vector{T2}} where T2<:Number)
        @boundscheck _chkbounds_Pl(pl, lmax, x)
        const U = promote_type(eltype(pl), eltype(x), Float64)

        @inbounds begin
            pl[:,1] .= 1.0
            lmax == 0 && return pl

            pl[:,2] .= x
            lmax == 1 && return pl

            for ℓ=2:lmax
                const onebyℓ = one(U)/ℓ
                pl[:,ℓ+1] .= (2.0.-onebyℓ).*x.*pl[:,ℓ] .- (1.0.-onebyℓ).*pl[:,ℓ-1]
            end
        end

        return pl
    end

    #########################################################################
    # ASSOCIATED LEGENDRE FUNCTIONS
    #########################################################################

    # Do expanded-precision calculations of several initalization constants
    # with the expectation that reducing the error on the first terms will
    # result in less error for all terms.
    const bigpi = convert(BigFloat, π);
    const sqrt1by2pi = convert(Float64, sqrt(big"0.50"/bigpi));
    const sqrt3by2pi = convert(Float64, sqrt(big"1.50"/bigpi));
    const sqrt3by4pi = convert(Float64, sqrt(big"0.75"/bigpi));

    doc"""
        immutable PlmPlan{T<:AbstractFloat}

    Structure containing precomputed coefficients for use with `Plm` and
    `Plm!` when computing the associated Legendre functions.
    """
    immutable PlmPlan{T<:AbstractFloat,L}
        # Pre-computed coefficients which are common to all calculates of P_l^m
        # no matter what the argument value is.
        #
        # Note that due to the form of the recursion relations available, we'll
        # never need the coefficients for any of (ℓ,m)∈{(0,0),(1,0),(1,1)}, nor
        # for any m==ℓ-1 and m==ℓ.
        α::Matrix{T}
        β::Matrix{T}
        # Additionally, there are m-only coefficients which are repeated:
        μ::Vector{T}
        ν::Vector{T}
        # To allow code reuse for multiple normalization conventions, we
        # store the coefficients for the ℓ=0 and ℓ=1 initial conditions.
        Plm00::T
        Plm10::T
        Plm11::T
    end

    doc"""
        PlmPlan{T<:AbstractFloat}(lmax::Integer)

    Returns a `PlmPlan` structure which has precomputed coefficients for
    associated Legendre functions up to degree `lmax` using the unit
    normalization convention.
    """
    function PlmPlan(T::Type, lmax::Integer)
        # Use doubles for internal computations unless BigFloat is used,
        # in which case we want to store extra precision.
        const U = promote_type(T, Float64)

        const Nm = lmax + 1
        α = Matrix{T}(Nm, Nm)
        β = Matrix{T}(Nm, Nm)
        μ = Vector{T}(Nm)
        ν = Vector{T}(Nm)

        # Loop over all ℓ > 1
        @inbounds for ℓ=2:lmax
            const ℓℓm1 = convert(U, 2*ℓ-1)
            const ℓm1 = convert(U, ℓ-1)

            μ[ℓ+1] = -ℓℓm1
            ν[ℓ+1] =  ℓℓm1

            for m=0:(ℓ-2)
                α[m+1,ℓ+1] = ℓℓm1 / (ℓ-m)
                β[m+1,ℓ+1] = (ℓm1+m) / ℓℓm1
            end
        end

        # Deal with the initial conditions:
        #
        # P_0^0 =  1            =  1
        # P_1^0 =  cos(θ)*P_0^0 =  cos(θ)
        # P_1^1 = -sin(θ)*P_0^0 = -sin(θ)
        #
        # The leading coefficients have already been calculated into
        # the constants during initialization.
        PlmPlan{T,lmax}(α, β, μ, ν, T(1), T(1), T(-1))
    end

    PlmPlan(L::Integer) = PlmPlan(Float64, L);

    doc"""
        PlmSphericalPlan(T::Type, lmax::Integer)

    Returns a `PlmPlan` structure which has precomputed coefficients for
    associated Legendre functions up to degree `lmax` such that the spherical
    harmonics $Y_l^m$ can be efficiently computed.
    """
    function PlmSphericalPlan(T::Type, lmax::Integer)
        # Use doubles for internal computations unless BigFloat is used,
        # in which case we want to store extra precision.
        const U = promote_type(T, Float64)

        const Nm = lmax + 1
        α = Matrix{T}(Nm, Nm)
        β = Matrix{T}(Nm, Nm)
        μ = Vector{T}(Nm)
        ν = Vector{T}(Nm)

        # Loop over all ℓ > 1
        @inbounds for ℓ=2:lmax
            const ℓd = convert(U, ℓ)
            const ℓsq = ℓd * ℓd
            const ℓm1sq = (ℓd-U(1))*(ℓd-U(1))

            μ[ℓ+1] = -sqrt(1 + U(0.5)/ℓ)
            ν[ℓ+1] =  sqrt(U(2)*ℓ + 1)

            for m=0:(ℓ-2)
                const msq = convert(U,m) * convert(U,m)
                α[m+1,ℓ+1] = sqrt( (U(4)*ℓsq - U(1)) / (ℓsq - msq) )
                β[m+1,ℓ+1] = sqrt( (ℓm1sq - msq) / (U(4)*ℓm1sq - U(1)) )
            end
        end

        # Deal with the initial conditions:
        #
        # P_0^0 = sqrt(2) * sqrt(1/4π)    =  sqrt(1/2π)
        # P_1^0 =  sqrt(3)*cos(θ)*P_0^0   =  sqrt(3/2π) * cos(θ)
        # P_1^1 = -sqrt(3/2)*sin(θ)*P_0^0 = -sqrt(3/4π) * sin(θ)
        #
        # The leading coefficients have already been calculated into
        # the constants during initialization.
        PlmPlan{T,lmax}(α, β, μ, ν, sqrt1by2pi, sqrt3by2pi, -sqrt3by4pi)
    end

    PlmSphericalPlan(lmax::Integer) = PlmSphericalPlan(Float64, lmax)

    doc"""
        Plm{L}(pp::PlmPlan{T,L} where T, x)

    Computes normalized associated Legendre function values at `x` for all
    $0 \le \ell \le L$, $0 \le m \le l$. `pp` is a design structure containing
    precomputed values for a given max `L` which are constant for all values of
    `x`.
    """
    function Plm end

    function Plm{L}(pp::PlmPlan{T1,L} where T1<:Number, x::Number)
        plm = Array{eltype(x)}(1, L+1, L+1)
        @inbounds Plm!(plm, pp, x)
        return reshape(plm, L+1, L+1)
    end
    function Plm{L}(pp::PlmPlan{T1,L} where T1<:Number, x::Vector{T2} where T2<:Number)
        plm = Array{eltype(x)}(length(x), L+1, L+1)
        @inbounds Plm!(plm, pp, x)
        return plm
    end

    @inline function _chkbounds_Plm(plm, lmax, x)
        lmax >= 0 || throw(DimensionMismatch("lmax must be non-negative"))
        size(plm) >= (length(x),lmax+1,lmax+1) || throw(DimensionMismatch())
    end

    doc"""
        Plm!{L}(plm, pp::PlmPlan{T,L} where T, x)

    Computes normalized associated Legendre function values at `x` for all
    $0 \le \ell \le L$, $0 \le m \le l$. The answer is stored in `plm`.
    `pp` is a design structure containing precomputed values for a
    given lmax `L` (which are constant for all values of `x`).
    """
    @propagate_inbounds function Plm!{T<:Number,L}(plm::Array{T1,3} where T1<:Number, pp::PlmPlan{T,L}, x::Union{T2,Vector{T2}} where T2<:Number)
        @boundscheck _chkbounds_Plm(plm, L, x)
        const U = promote_type(eltype(plm), T, eltype(x), Float64)

        # Given that x = cos(θ), we can get sin(θ) by trig identities:
        const y = sqrt.(one(U) .- x.*x)

        # To start the recusion, we need ℓ=0 and ℓ=1 filled in. We know these
        # analytically:
        @inbounds plm[:,1,1] .= pp.Plm00
        L == 0 && return plm

        @inbounds plm[:,1,2] .= pp.Plm10 .* x
        @inbounds plm[:,2,2] .= pp.Plm11 .* y
        L == 1 && return plm

        @inbounds for ℓ=2:L
            # Do most of the terms using the recursion relation
            #
            #   P_ℓ^m = α_ℓ^m * (x*P_(ℓ-1)^m - β_ℓ^m * P_(ℓ-2)^m)
            #
            # where α and β are defined by the design structure.
            #
            # Don't forget that indexing needs + 1 since ℓ=0/m=0 is element 1!
            for m=0:(ℓ-2)
                plm[:,m+1,ℓ+1] .= pp.α[m+1,ℓ+1] .* (x.*plm[:,m+1,ℓ] .-
                        pp.β[m+1,ℓ+1].*plm[:,m+1,ℓ-1])
            end

            # Then for the last two terms, we use new recursion relations
            # where again the leading coefficients have been pre-computed
            # and stored in the design structure:
            #
            #   P_ℓ^(ℓ-1) = ν_ℓ * x * P_(ℓ-1)^(ℓ-1)
            #
            plm[:,ℓ,ℓ+1] .= pp.ν[ℓ+1] .* x .* plm[:,ℓ,ℓ]
            #
            # and
            #
            #   P_ℓ^ℓ = μ_ℓ * y * P_(ℓ-1)^(ℓ-1)
            #
            plm[:,ℓ+1,ℓ+1] .= pp.μ[ℓ+1] .* y .* plm[:,ℓ,ℓ]
        end

        return plm
    end
end

