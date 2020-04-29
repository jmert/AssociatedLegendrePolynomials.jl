using Random
import ..Norms

@testset "Equality of legendre[!]" begin
    LMAX = 10
    x = 0.5
    λ₂ = fill(0.0, LMAX+1)
    Λ₂ = fill(0.0, LMAX+1, LMAX+1)

    leg! = LegendreSphereCoeff{Float64}(LMAX)
    # Fill lower triangular matrix by hand
    Λ₁ = [m > l ? 0.0 : λlm(l, m, x) for l in 0:LMAX, m in 0:LMAX]

    # Test full matrix
    @test Λ₁ == leg!(Λ₂, LMAX, LMAX, x)
    @test Λ₁ == λlm!(Λ₂, LMAX, LMAX, x)

    # Test single columns
    fill!(λ₂, 0.0)
    @test @view(Λ₁[:,2+1]) == leg!(λ₂, LMAX, 2, x)
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

@testset "Return inference ($(typeof(N)) with argument type $T)" for N in Norms, T in NumTypes
    import ..LMAX
    Λ₀ = fill(0.0)
    Λ₁ = zeros(LMAX+1)
    Λ₂ = zeros(LMAX+1, LMAX+1)
    # 7 cases in recursion to handle:
    #   1. Initial condition
    #   2. Boost along diagonal from (even,even) -> (odd,odd)
    #      2a. Boost in ℓ from (m,m) -> (m+1, m)
    #      2b. Boost again in ℓ from (m+1,m) -> (m+2, m)
    #   3. Repeat (2) but for step along diagonal from (odd,odd) -> (even,even)
    cases = ((0,0), (1,1), (2,1), (3,1), (2,2), (3,2), (4,2))
    @testset "(ℓ,m) == ($ℓ, $m)" for (ℓ, m) in cases
        # Return type matches that of the argument, even though internal calculation
        # will promote.
        @test @inferred(legendre(N, ℓ, m, one(T))) isa T
        # Also must handle single value, all ℓ for fixed m, and all (ℓ,m) up to
        # (LMAX,LMAX)
        @test @inferred(legendre!(N, Λ₀, LMAX, LMAX, one(T))) isa typeof(Λ₀)
        @test @inferred(legendre!(N, Λ₁, LMAX, LMAX, one(T))) isa typeof(Λ₁)
        @test @inferred(legendre!(N, Λ₂, LMAX, LMAX, one(T))) isa typeof(Λ₂)
    end
end

@testset "Internal type promotion" begin
    LMAX = 5
    # The implementation should promote all internal calculations, so given arguments
    # 0.5f0 (32-bit) and 0.5e0 (64-bit) which can both be exactly represented, calculate
    # that the former is calculated at 64-bit precision when the 64-bit output encourages
    # it.
    @test λlm!(fill(0.0), LMAX, LMAX, 0.5f0) == λlm!(fill(0.0), LMAX, LMAX, 0.5e0)
    # Do the same for a normalization which contains its own inherent element type.
    # Now take advantage of the fact that 0.2 is not exactly representable and have
    # different approximations at different accuracy to induce numerical differences
    # based on whether the values have been promoted or not.
    bigΛ! = LegendreSphereCoeff{BigFloat}(LMAX)
    λ1 = bigΛ!(fill(0f0), LMAX, LMAX, 0.2f0)
    λ2 = bigΛ!(fill(big(0f0)), LMAX, LMAX, big(0.2f0))
    @test λ1[] == Float32(λ2[])
end

@testset "Broadcasting arguments" begin
    LMAX = 5
    ctab = LegendreSphereCoeff{Float64}(LMAX)
    x = 0.5

    # Single (l,m) for single z
    @test Plm.(LMAX, LMAX, x) == Plm(LMAX, LMAX, x)
    @test λlm.(LMAX, LMAX, x)  == λlm(LMAX, LMAX, x)
    @test ctab.(LMAX, LMAX, x) == λlm(LMAX, LMAX, x)
    @test legendre.(LegendreSphereNorm(), LMAX, LMAX, x) == λlm(LMAX, LMAX, x)

    z = range(-1, 1, length=10)
    Λ = zeros(length(z), LMAX + 1, LMAX + 1)

    # Single (l,m) over multiple z
    legendre!(LegendreUnitNorm(), Λ, LMAX, LMAX, z)
    @test Plm.(LMAX, LMAX, z) == Λ[:,LMAX+1,LMAX+1]
    legendre!(LegendreSphereNorm(), Λ, LMAX, LMAX, z)
    @test λlm.(LMAX, LMAX, z)  == Λ[:,LMAX+1,LMAX+1]
    @test ctab.(LMAX, LMAX, z) == Λ[:,LMAX+1,LMAX+1]
    @test legendre.(LegendreSphereNorm(), LMAX, LMAX, z) == Λ[:,LMAX+1,LMAX+1]

    # All l for fixed m over multiple z
    legendre!(LegendreUnitNorm(), Λ, LMAX, LMAX, z)
    @test Plm.(0:LMAX, LMAX, z) == Λ[:,:,LMAX+1]
    legendre!(LegendreSphereNorm(), Λ, LMAX, LMAX, z)
    @test λlm.(0:LMAX, LMAX, z)  == Λ[:,:,LMAX+1]
    @test ctab.(0:LMAX, LMAX, z) == Λ[:,:,LMAX+1]
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
end
