using LinearAlgebra, Random

import ..TestSuite
TestSuite.runtests(LegendreUnitNorm())
TestSuite.runtests(LegendreSphereNorm())
TestSuite.runtests(LegendreUnitCoeff{Float64}(5))

@testset "Coefficient table conversion" begin
    dtab = LegendreSphereCoeff{Float64}(10)
    ftab = convert(LegendreSphereCoeff{Float32}, dtab)

    @test ftab isa LegendreSphereCoeff{Float32}
    @test ftab.α == Float32.(dtab.α)
    @test_throws MethodError convert(LegendreUnitCoeff{Float64}, dtab)
end

@testset "Setting mmax" begin
    LMAX = 10
    MMAX = 2
    impltab = LegendreUnitCoeff{Float64}(LMAX)
    fulltab = LegendreUnitCoeff{Float64}(LMAX, LMAX)
    mmaxtab = LegendreUnitCoeff{Float64}(LMAX, MMAX)

    # Equality of implicit and explicit [maximal] mmax
    @test impltab.α == fulltab.α
    @test impltab.α == fulltab.α
    @test impltab.β == fulltab.β
    @test impltab.μ == fulltab.μ
    @test impltab.ν == fulltab.ν

    # Equality of coefficients up to mmax for mmax limited
    @test mmaxtab.α == @view fulltab.α[:, 1:(MMAX+1)]
    @test mmaxtab.β == @view fulltab.β[:, 1:(MMAX+1)]
    @test mmaxtab.μ == @view fulltab.μ[1:(MMAX+1)]
    @test mmaxtab.ν == @view fulltab.ν[1:(MMAX+1)]
end

@testset "Domain checking" begin
    LMAX = 10
    MMAX = 2
    ctab = LegendreUnitCoeff{Float64}(LMAX)
    mtab = LegendreUnitCoeff{Float64}(LMAX, MMAX)
    λ = Vector{Float64}(undef, LMAX+1)
    Λ = Matrix{Float64}(undef, LMAX+1, LMAX+1)

    # Mathematical domain errors:
    @test_throws DomainError Pl(-1, 0.5)
    @test_throws DomainError Plm(-1, 0, 0.5)
    @test_throws DomainError Plm(LMAX, -1, 0.5)
    @test_throws DomainError Plm(LMAX, LMAX+1, 0.5)
    @test_throws DomainError Pl!(λ, -1, 0.5)
    @test_throws DomainError Plm!(Λ, -1, 0, 0.5)
    @test_throws DomainError Plm!(Λ, LMAX, -1, 0.5)
    @test_throws DomainError Plm!(Λ, LMAX, LMAX+1, 0.5)
    @test_throws DomainError legendre(ctab, -1, 0, 0.5)
    @test_throws DomainError legendre(ctab, LMAX, -1, 0.5)
    @test_throws DomainError legendre(ctab, LMAX, LMAX+1, 0.5)
    @test_throws DomainError legendre!(ctab, λ, -1, LMAX, 0.5)
    @test_throws DomainError legendre!(ctab, Λ, -1, LMAX, 0.5)
    @test_throws DomainError legendre!(ctab, Λ, LMAX, -1, 0.5)
    @test_throws DomainError legendre!(ctab, Λ, LMAX, LMAX+1, 0.5)

    # Bounds error for precomputed coefficient tables
    @test_throws BoundsError legendre(ctab, LMAX+1, 0, 0.5)
    @test_throws BoundsError legendre(mtab, LMAX, MMAX+1, 0.5)
    @test_throws BoundsError legendre!(mtab, λ, LMAX, MMAX+1, 0.5)
    @test_throws BoundsError legendre!(mtab, Λ, LMAX, MMAX+1, 0.5)
end

@testset "Output array bounds checking" begin
    LMAX = 2
    ctab = LegendreUnitCoeff{Float64}(LMAX)
    λ  = Vector{Float64}(undef, LMAX)
    Λ₁ = Matrix{Float64}(undef, LMAX, LMAX+1)
    Λ₂ = Matrix{Float64}(undef, LMAX+1, LMAX)

    # Bounds error for filling vector or matrix
    @test_throws DimensionMismatch Pl!(λ, LMAX, 0.5)
    @test_throws DimensionMismatch Plm!(λ, LMAX, 0, 0.5)
    @test_throws DimensionMismatch Plm!(Λ₁, LMAX, LMAX, 0.5)
    @test_throws DimensionMismatch Plm!(Λ₂, LMAX, LMAX, 0.5)
    @test_throws DimensionMismatch legendre!(ctab, λ, LMAX, 0, 0.5)
    @test_throws DimensionMismatch legendre!(ctab, Λ₁, LMAX, LMAX, 0.5)
    @test_throws DimensionMismatch legendre!(ctab, Λ₂, LMAX, LMAX, 0.5)

    # Compatibility of output storage
    Λ₃ = Array{Float64}(undef, 2, LMAX+1, LMAX+1)
    Λ₄ = Array{Float64}(undef, 2, 2, LMAX+1, LMAX+1)
    Λ₅ = Array{Float64}(undef, 2, 2, 2, LMAX+1, LMAX+1)
    # Insufficient dimensions:
    @test_throws DimensionMismatch Plm!(Λ₁, LMAX, LMAX, 0.0)
    @test_throws DimensionMismatch Plm!(Λ₂, LMAX, LMAX, zeros(2))
    @test_throws DimensionMismatch Plm!(Λ₃, LMAX, LMAX, zeros(2,2))
    # Too many dimensions:
    @test_throws DimensionMismatch Plm!(Λ₃, LMAX, LMAX, 0.0)
    @test_throws DimensionMismatch Plm!(Λ₄, LMAX, LMAX, zeros(2))
    @test_throws DimensionMismatch Plm!(Λ₅, LMAX, LMAX, zeros(2,2))
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
    @test all(leg(Λ₁, 0.5) .== legendre!(leg, Λ₂, LMAX, LMAX, 0.5))
end

@testset "Mixed types" begin
    LMAX = 10
    legb! = LegendreUnitCoeff{BigFloat}(LMAX)
    legd! = LegendreUnitCoeff{Float64}(LMAX)
    λb = Vector{BigFloat}(undef, LMAX+1)
    λd = Vector{Float64}(undef, LMAX+1)
    xb = big"0.5"
    xd = 5e-1

    # Same table and array, mixed value
    @test @inferred(legb!(λb, 2, xd)) isa typeof(λb)
    @test @inferred(legd!(λd, 2, xb)) isa typeof(λd)

    # Same table and value, mixed array
    @test @inferred(legb!(λd, 2, xb)) isa typeof(λd)
    @test @inferred(legd!(λb, 2, xd)) isa typeof(λb)

    # Same array and value, mixed table
    @test @inferred(legb!(λd, 2, xd)) isa typeof(λd)
    @test @inferred(legd!(λb, 2, xb)) isa typeof(λb)

    # All three mixed
    @test @inferred(legb!(λd, 2, Float32(xd))) isa typeof(λd)
end

@testset "Equality of legendre[!]" begin
    LMAX = 10
    x = 0.5
    λ₂ = fill(0.0, LMAX+1)
    Λ₂ = fill(0.0, LMAX+1, LMAX+1)

    # P_l(x) is implemented with its own fast-path for m == 0
    λ₁ = [Pl(l, x) for l in 0:LMAX]
    @test λ₁ == Pl!(λ₂, LMAX, x)

    # Use λlm instead of Plm to provide coverage for the convenience wrapper functions,
    # too
    leg! = LegendreSphereCoeff{Float64}(LMAX)
    # Fill lower triangular matrix by hand
    Λ₁ = [m > l ? 0.0 : λlm(l, m, x) for l in 0:LMAX, m in 0:LMAX]

    # Test full matrix
    @test Λ₁ == leg!(Λ₂, x)
    @test Λ₁ == λlm!(Λ₂, LMAX, LMAX, x)

    # Test single columns
    fill!(λ₂, 0.0)
    @test @view(Λ₁[:,2+1]) == leg!(λ₂, 2, x)
    fill!(λ₂, 0.0)
    @test @view(Λ₁[:,2+1]) == λlm!(λ₂, LMAX, 2, x)
end

# Test failed with world age problems when generated functions were incorrectly used;
# numerical derivative via newer world's Dual numbers tests for this error.
@testset "Derivatives via dual numbers" begin
    using ForwardDiff: derivative
    LMAX = 5
    ctab = LegendreUnitCoeff{Float64}(LMAX)

    # Numerical derivative for Plm via recurrence relations:
    function num_deriv1(l, m, x)
        d = -(l+1) * x * Plm(l, m, x) + (l - m + 1) * Plm(l+1, m, x)
        return d / (x^2 - 1)
    end
    dual_deriv1(l, m, x) = derivative(z -> Plm(l, m, z), x)
    dual_deriv1_tab(l, m, x) = derivative(z -> ctab(l, m, z), x)

    for l in 0:LMAX, m in 0:l
        x = 2rand() - 1
        @test num_deriv1(l, m, x) ≈ dual_deriv1(l, m, x)
        @test num_deriv1(l, m, x) ≈ dual_deriv1_tab(l, m, x)
    end
end

@testset "Broadcasting arguments" begin
    LMAX = 5
    ctab = LegendreSphereCoeff{Float64}(LMAX)
    x = 0.5

    # Single (l,m) for single z
    @test Pl.(LMAX, x)        == Pl(LMAX, x)
    @test Plm.(LMAX, LMAX, x) == Plm(LMAX, LMAX, x)
    @test λlm.(LMAX, LMAX, x)  == λlm(LMAX, LMAX, x)
    @test ctab.(LMAX, x)       == λlm(LMAX, 0, x)
    @test ctab.(LMAX, LMAX, x) == λlm(LMAX, LMAX, x)
    @test legendre.(LegendreSphereNorm(), LMAX, x)       == λlm(LMAX, 0, x)
    @test legendre.(LegendreSphereNorm(), LMAX, LMAX, x) == λlm(LMAX, LMAX, x)

    z = range(-1, 1, length=10)
    Λ = zeros(length(z), LMAX + 1, LMAX + 1)

    # Single (l,m) over multiple z
    legendre!(LegendreUnitNorm(), Λ, LMAX, LMAX, z)
    @test Pl.(LMAX, z)        == Λ[:,LMAX+1,1]
    @test Plm.(LMAX, LMAX, z) == Λ[:,LMAX+1,LMAX+1]
    legendre!(LegendreSphereNorm(), Λ, LMAX, LMAX, z)
    @test λlm.(LMAX, LMAX, z)  == Λ[:,LMAX+1,LMAX+1]
    @test ctab.(LMAX, z)       == Λ[:,LMAX+1,1]
    @test ctab.(LMAX, LMAX, z) == Λ[:,LMAX+1,LMAX+1]
    @test legendre.(LegendreSphereNorm(), LMAX, z)       == Λ[:,LMAX+1,1]
    @test legendre.(LegendreSphereNorm(), LMAX, LMAX, z) == Λ[:,LMAX+1,LMAX+1]

    # All l for fixed m over multiple z
    legendre!(LegendreUnitNorm(), Λ, LMAX, LMAX, z)
    @test Pl.(0:LMAX, z)        == Λ[:,:,1]
    @test Plm.(0:LMAX, LMAX, z) == Λ[:,:,LMAX+1]
    legendre!(LegendreSphereNorm(), Λ, LMAX, LMAX, z)
    @test λlm.(0:LMAX, LMAX, z)  == Λ[:,:,LMAX+1]
    @test ctab.(0:LMAX, z)       == Λ[:,:,1]
    @test ctab.(0:LMAX, LMAX, z) == Λ[:,:,LMAX+1]
    @test legendre.(LegendreSphereNorm(), 0:LMAX, z)       == Λ[:,:,1]
    @test legendre.(LegendreSphereNorm(), 0:LMAX, LMAX, z) == Λ[:,:,LMAX+1]

    # All l and m over multiple z
    legendre!(LegendreUnitNorm(), Λ, LMAX, LMAX, z)
    @test Plm.(0:LMAX, 0:LMAX, z) == Λ
    legendre!(LegendreSphereNorm(), Λ, LMAX, LMAX, z)
    @test λlm.(0:LMAX, 0:LMAX, z)  == Λ
    @test ctab.(0:LMAX, 0:LMAX, z) == Λ
    @test legendre.(LegendreSphereNorm(), 0:LMAX, 0:LMAX, z) == Λ

    # Shape-preservation of multi-dimensional arguments
    sz = (5, 10)
    z = collect(reshape(range(-1.0, 1.0, length=prod(sz)), sz))
    @test size(Plm.(LMAX, 0, z)) == sz
    @test size(Plm.(0:LMAX, 0, z)) == (sz..., LMAX+1)
    @test size(Plm.(0:LMAX, 0:LMAX, z)) == (sz..., LMAX+1, LMAX+1)

    # Throws on invalid l ranges
    @test_throws ArgumentError Pl.(1:LMAX, z)
    @test_throws ArgumentError Plm.(1:LMAX, 0, z)
    @test_throws ArgumentError λlm.(1:LMAX, 0, z)
    @test_throws ArgumentError ctab.(1:LMAX, 0, z)
    @test_throws ArgumentError legendre.(ctab, 1:LMAX, 0, z)
    # Throws on invalid m ranges
    @test_throws ArgumentError Plm.(0:LMAX, 1:LMAX, z)
    @test_throws ArgumentError λlm.(0:LMAX, 1:LMAX, z)
    @test_throws ArgumentError ctab.(0:LMAX, 1:LMAX, z)
    @test_throws ArgumentError legendre.(ctab, 0:LMAX, 1:LMAX, z)
end

@testset "Numerical stability (issue #11)" begin
    x = sind(-57.5)
    lmax = 3 * 720
    Λ = λlm.(0:lmax, 0:lmax, x)
    # The amplitude bound of 1.5 is a very rough limit --- simply checking for
    # unbounded growth.
    @test all(abs.(Λ[end,:]) .< 1.5)
end
