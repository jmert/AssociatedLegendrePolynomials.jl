using LinearAlgebra

################################
# ASSOCIATED LEGENDRE FUNCTIONS
################################

# P_0^0 is constant. Verify the output is invariant for several inputs.
@testset "Constant P_0^0 ($T)" for T in NumTypes
    @test @inferred(legendre(LegendreUnitNorm(), 0, 0, T(0.1))) isa T
    @test all(legendre(LegendreUnitNorm(), 0, 0, z) == one(T)
              for z in range(-one(T), one(T), length=10))
end

# Match P_ℓ^{m=0} terms for ℓ=1...9 where the analytical expressions
# have been taken from a table of equations, assuming x = cos(π/4).
@testset "Analytic checks for P_ℓ^0(cos(π/4)) ($T)" for T in NumTypes
    LMAX = 9
    plm = Matrix{T}(undef, LMAX+1, LMAX+1)
    legendre!(LegendreUnitNorm(), plm, LMAX, LMAX, sqrt(T(2))/2)
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

##############################################################
# SPHERICAL HARMONIC NORMALIZED ASSOCIATED LEGENDRE FUNCTIONS
##############################################################

# P_0^0 is constant. Verify the output is invariant for several inputs.
@testset "Constant P_0^0 ($T)" for T in NumTypes
    @test @inferred(legendre(LegendreSphereNorm(), 0, 0, T(0.1))) isa T
    @test all(legendre(LegendreSphereNorm(), 0, 0, z) == sqrt(inv(4T(π)))
              for z in range(-one(T), one(T), length=10))
end

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
        ratio(x) = sqrt(T(2x + 1) / T(2x))
        factorials = mapreduce(ratio, *, 1:ℓ, init=one(T))
        sign = (ℓ % 2 == 0) ? 1 : -1
        return sign * sqrt(inv(4T(π))) * factorials
    end

    LMAX = 99
    Λ = fill(0.0, LMAX+1, LMAX+1)
    legendre!(LegendreSphereNorm(), Λ, LMAX, LMAX, cosd(45.0))

    @test Λ[3,3]   ≈ SphericalPll_coeff(Float64,2)  / (2^1)
    @test Λ[4,4]   ≈ SphericalPll_coeff(Float64,3)  / (2^1 * sqrt(2))
    @test Λ[21,21] ≈ SphericalPll_coeff(Float64,20) / (2^10)
    @test Λ[22,22] ≈ SphericalPll_coeff(Float64,21) / (2^10 * sqrt(2))
    @test Λ[99,99] ≈ SphericalPll_coeff(Float64,98) / (2^49)
    @test Λ[100,100] ≈ SphericalPll_coeff(Float64,99) / (2^49 * sqrt(2))
end

# The normal P_ℓ^m should be equal to the spherical harmonic normalized
# λ_ℓ^m if we manually normalize them.
@testset "Equality of N_ℓ^m × P_ℓ^m and λ_ℓ^m ($T)" for T in NumTypes
    LMAX = 5
    atol = max(eps(T(100)), eps(100.0)) # maximum(plm_norm) is of order 10²
    plm_norm = zeros(T, LMAX+1, LMAX+1)
    plm_sphr = zeros(T, LMAX+1, LMAX+1)

    lmat = tril(repeat(collect(0:LMAX), 1, LMAX+1))
    mmat = tril(repeat(collect(0:LMAX)', LMAX+1, 1))
    nlm = tril(Nlm.(T, lmat, mmat))
    for ii in 1:10
        x = 2*T(rand()) - 1
        legendre!(LegendreUnitNorm(), plm_norm, LMAX, LMAX, x)
        legendre!(LegendreSphereNorm(), plm_sphr, LMAX, LMAX, x)
        @test all(isapprox.(nlm.*plm_norm, plm_sphr, atol=atol))
    end
end

################################
# Complex domain
################################

# Trivial case: real-only complex arguments should be identical to the real case
@testset "Complex arguments, real axis ($T)" for T in NumTypes
    LMAX = 10
    x = collect(range(-one(T), one(T), length=100))
    z = complex(x)
    @test Plm.(0:LMAX, 0:LMAX, x) == real.(Plm.(0:LMAX, 0:LMAX, z))
    @test λlm.(0:LMAX, 0:LMAX, x) == real.(λlm.(0:LMAX, 0:LMAX, z))
end

################################
# Other analytic checks / issues
################################

@testset "Numerical stability/accuracy" begin
    # issue CMB.jl#11
    #   Overall accuracy --- fixed by adding FMA y=(1-x^2) and pairwise diagonal iteration
    x = sind(-57.5)
    lmax = 3 * 720
    Λ = λlm.(0:lmax, 0:lmax, x)
    # The amplitude bound of 1.5 is a very rough limit --- simply checking for
    # unbounded growth.
    @test all(abs.(Λ[end,:]) .< 1.5)

    # PR#16
    #   Increase accuracy of μ coefficient for spherical norm.
    #   Also test ν, and both already equivalent for unit norm.
    @testset "coeffs μ, ν ($norm)" for norm in (LegendreSphereNorm(), LegendreUnitNorm())
        lrng = 1:30_000
        μ(T, l) = Legendre.coeff_μ(norm, T, l)
        ν(T, l) = Legendre.coeff_ν(norm, T, l)
        @test μ.(Float64, lrng) == Float64.(μ.(BigFloat, lrng))
        @test ν.(Float64, lrng) == Float64.(ν.(BigFloat, lrng))
    end
end
