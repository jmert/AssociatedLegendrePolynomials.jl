module Legendre
    using Compat.Test
    using ..CMBTests: NumTypes
    using CMB.Legendre

    @testset "Domain and bounds checking" begin
        LMAX = 10
        ctab = LegendreUnitCoeff{Float64}(LMAX)

        # Mathematical domain errors:
        @test_throws DomainError Plm(-1, 0, 0.5)
        @test_throws DomainError Plm(LMAX, -1, 0.5)
        @test_throws DomainError Plm(LMAX, LMAX+1, 0.5)
        @test_throws DomainError legendre(ctab, -1, 0, 0.5)
        @test_throws DomainError legendre(ctab, LMAX, -1, 0.5)
        @test_throws DomainError legendre(ctab, LMAX, LMAX+1, 0.5)

        # Bounds error for precomputed coefficient tables
        @test_throws BoundsError legendre(ctab, LMAX+1, 0, 0.5)

        # Bounds error for filling vector or matrix
        λ = Vector{Float64}(LMAX)
        Λ₁ = Matrix{Float64}(LMAX, LMAX+1)
        Λ₂ = Matrix{Float64}(LMAX+1, LMAX)
        @test_throws DimensionMismatch legendre!(ctab, λ, LMAX, 0, 0.5)
        @test_throws DimensionMismatch legendre!(ctab, Λ₁, LMAX, 0.5)
        @test_throws DimensionMismatch legendre!(ctab, Λ₂, LMAX, 0.5)
    end

    @testset "Functor interface" begin
        LMAX = 10
        leg = LegendreSphereCoeff{Float64}(LMAX)
        λ₁ = fill(0.0, LMAX+1)
        λ₂ = fill(0.0, LMAX+1)
        Λ₁ = fill(0.0, LMAX+1, LMAX+1)
        Λ₂ = fill(0.0, LMAX+1, LMAX+1)
        @test leg(1, 1, 0.5) == legendre(leg, 1, 1, 0.5)
        @test leg(1, 0.5) == legendre(leg, 1, 0.5)
        @test all(leg(λ₁, 2, 0.5) .== legendre!(leg, λ₂, LMAX, 2, 0.5))
        @test all(leg(Λ₁, 0.5) .== legendre!(leg, Λ₂, LMAX, 0.5))
    end

    #######################
    # LEGENDRE POLYNOMIALS
    #######################

    srand(2222)

    # P_0 is a constant. Verify the output is invariant.
    @testset "Constant P_0 ($T)" for T in NumTypes
        @test @inferred(Pl(0, T(0.1))) isa T
        pl = Pl.(0, 2 .* T.(rand(10)) .- 1)
        @test all(pl[:,1] .== T(1.0))
    end

    ################################
    # ASSOCIATED LEGENDRE FUNCTIONS
    ################################

    # P_0^0 is constant. Verify the output is invariant for several inputs.
    @testset "Constant P_0^0 ($T)" for T in NumTypes
        @test @inferred(Plm(0, 0, T(0.1))) isa T
        plm = Plm.(0, 0, 2 .* T.(rand(10)) .- 1)
        @test all(plm .== T(1.0))
    end

    # Match P_ℓ^{m=0} terms for ℓ=1...9 where the analytical expressions
    # have been taken from a table of equations, assuming x = cos(π/4).
    @testset "Analytic checks for P_ℓ^0(cos(π/4)) ($T)" for T in NumTypes
        LMAX = 9
        ctab = LegendreUnitCoeff{T}(LMAX)
        plm = Matrix{T}(LMAX+1, LMAX+1)
        legendre!(ctab, plm, LMAX, sqrt(T(2))/2)
        @test plm[2, 1] ≈  sqrt(T(2))/2
        @test plm[3, 1] ≈  T(1)/4
        @test plm[4, 1] ≈ -sqrt(T(2))/8
        @test plm[5, 1] ≈ -T(13)/32
        @test plm[6, 1] ≈ -17sqrt(T(2))/64
        @test plm[7, 1] ≈ -T(19)/128
        @test plm[8, 1] ≈  23sqrt(T(2))/256
        @test plm[9, 1] ≈  T(611)/2048
        @test plm[10,1] ≈  827sqrt(T(2))/4096
    end

    # The associated Legendre functions P_ℓ^m should equate to the Legendre
    # polynomials P_ℓ for m=0, so check this is true.
    @testset "Equality of P_ℓ and P_ℓ^0 ($T)" for T in NumTypes
        LMAX = 10
        ctab = LegendreUnitCoeff{T}(LMAX)
        pl  = Vector{T}(LMAX+1)
        plm = Matrix{T}(LMAX+1, LMAX+1)
        for ii in 1:10
            x = 2T(rand()) - 1
            Pl!(pl, LMAX, x)
            legendre!(ctab, plm, LMAX, x)
            @test all(pl .≈ plm[:,1])
        end
    end

    ##############################################################
    # SPHERICAL HARMONIC NORMALIZED ASSOCIATED LEGENDRE FUNCTIONS
    ##############################################################

    # The initialization condition
    #
    #   P_ℓ^ℓ(x) = (-1)^ℓ (2ℓ-1)!! (1-x^2)^(ℓ/2)
    #
    # can be used with extended-precision computations for special values
    # of ℓ and m. To account for the spherical harmonic normalization,
    # the formula can be simplified to
    #
    #   λ_ℓ^ℓ(x) = (-1)^ℓ 1/sqrt(4π) sqrt((2ℓ+1)!!/(2ℓ)!!) (1-x^2)^(ℓ/2)
    #
    # We'll compute the coefficient separately since it's "easy", and then
    # we can use convenient angles that also expand to an exact analytic
    # value for (1-x^2)^(ℓ/2).
    #
    # Assuming x = cos(θ), the (1-x^2)^(ℓ/2) term simplifies to an
    # expression of \sin(θ):
    #
    #   (1-x^2)^(ℓ/2) = (sin^2 θ)^(ℓ/2) = (sin θ)^ℓ
    #
    # For θ = π/4:
    #
    #   * even ℓ => 1 / 2^(ℓ/2)
    #   * odd ℓ => 1 / (2^(k/2) sqrt(2)) where k=floor(ℓ/2).
    @testset "Analytic checks for λ_ℓ^m(cos(π/4))" begin
        function SphericalPll_coeff(T, ℓ)
            factorials = mapreduce(x->sqrt((2x+1)/(2x)), *, 1.0, 1:ℓ)
            sign = (ℓ % 2 == 0) ? 1 : -1
            return sign * sqrt(1/4π) * factorials
        end

        LMAX = 99
        ctab = LegendreSphereCoeff{Float64}(LMAX)
        Λ = Matrix{Float64}(LMAX+1, LMAX+1)
        Λ = legendre!(ctab, Λ, LMAX, cosd(45.0))

        @test Λ[3,3]   ≈ SphericalPll_coeff(Float64,2)  / (2^1)
        @test Λ[4,4]   ≈ SphericalPll_coeff(Float64,3)  / (2^1 * sqrt(2))
        @test Λ[21,21] ≈ SphericalPll_coeff(Float64,20) / (2^10)
        @test Λ[22,22] ≈ SphericalPll_coeff(Float64,21) / (2^10 * sqrt(2))
        @test Λ[99,99] ≈ SphericalPll_coeff(Float64,98) / (2^49)
        @test Λ[100,100] ≈ SphericalPll_coeff(Float64,99) / (2^49 * sqrt(2))
    end

    # The normal P_ℓ^m should be equal to the spherical harmonic normalized
    # λ_ℓ^m if we manually normalize them.
    @testset "Equality of N_ℓ^m*P_ℓ^m and λ_ℓ^m ($T)" for T in NumTypes
        LMAX = 5
        atol = max(eps(T(100)), eps(100.0)) # maximum(plm_norm) is of order 10²
        ctab_norm = LegendreUnitCoeff{T}(LMAX)
        ctab_sphr = LegendreSphereCoeff{T}(LMAX)
        plm_norm = zeros(T, LMAX+1, LMAX+1)
        plm_sphr = zeros(T, LMAX+1, LMAX+1)

        lmat = tril(repmat(collect(0:LMAX), 1, LMAX+1))
        mmat = tril(repmat(collect(0:LMAX)', LMAX+1, 1))
        nlm = tril(Nlm.(T, lmat, mmat))
        for ii in 1:10
            x = 2*T(rand()) - 1
            legendre!(ctab_norm, plm_norm, LMAX, x)
            legendre!(ctab_sphr, plm_sphr, LMAX, x)
            @test all(isapprox.(nlm.*plm_norm, plm_sphr, atol=atol))
        end
    end
end

