using Random
using LinearAlgebra: tril, triu, diagind
import ..LMAX, ..Norms

@testset "Equality of legendre[!] ($(typeof(N)) with argument type $T)" for N in Norms, T in NumTypes
    x = T(0.5)
    # On the first pass, use out-of-place interface
    #
    # Fill the full matrix of values just once.
    λ2 = legendre(N, 0:LMAX, 0:LMAX, x)
    for m in 0:LMAX
        # Check that vector filling function returns same answers
        λ1 = legendre(N, 0:LMAX, m, x)
        @test @view(λ1[m+1:end]) == @view(λ2[m+1:end,m+1])
        for ℓ in m:LMAX
            # Then check that scalar output matches as well.
            @test legendre(N, ℓ, m, x) == λ1[ℓ+1]
        end
    end

    # Now on the second pass, use the in-place interface
    λ2 = fill(T(NaN), LMAX+1, LMAX+1)
    λ1 = fill(T(NaN), LMAX+1)
    λ0 = fill(T(NaN))

    legendre!(N, λ2, LMAX, LMAX, x)
    # Test for equality in legendre! and legendre; wrap with tril to exlude NaNs used later
    @test tril(λ2) == legendre(N, 0:LMAX, 0:LMAX, x)
    for m in 0:LMAX
        # Check that vector filling function returns same answers
        legendre!(N, λ1, LMAX, m, x)
        @test @view(λ1[m+1:end]) == @view(λ2[m+1:end,m+1])
        for ℓ in m:LMAX
            # Then check that scalar output matches as well.
            legendre!(N, λ0, ℓ, m, x)
            @test λ0[] == λ1[ℓ+1]
        end
    end

    # As an implementation detail, the outputs arrays are not cleared outside of the valid
    # domain when operated up by the in-place functions, so initial or previous values
    # should remain. Check these invariants as a way to notice an implementation change.
    #
    # For the matrix output, the upper triangle (excluding diagonal) should be all NaN.
    @test isequal(triu(λ2, 1), triu(fill(T(NaN), LMAX+1, LMAX+1), 1))
    # The vector output should accumulate the (m,m) term at each step, i.e. the diagonal
    # of the matrix output.
    @test λ1 == λ2[diagind(λ2)]
end

# Test failed with world age problems when generated functions were incorrectly used;
# numerical derivative via newer world's Dual numbers tests for this error.
@testset "Derivatives via dual numbers" begin
    using ForwardDiff: derivative
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
    Λ₀ = fill(0.0)
    Λ₁ = zeros(LMAX+1)
    Λ₂ = zeros(LMAX+1, LMAX+1)
    # 4 cases in recursion to handle:
    #   1. Initial condition
    #   2. Boost along diagonal
    #   3. Boost in ℓ from (m,m) -> (m+1, m)
    #   4. Boost again in ℓ from (m+1,m) -> (m+2, m) [and further]
    cases = ((0,0), (1,1), (2,1), (3,1))
    @testset "(ℓ,m) == ($ℓ, $m)" for (ℓ, m) in cases
        # Return type matches that of the argument for allocating interface, even though
        # internal calculation will promote.
        #
        # single (ℓ,m)
        @test @inferred(legendre(N, ℓ, m, one(T))) isa T
        @test @inferred(legendre!(N, Λ₀, LMAX, LMAX, one(T))) isa typeof(Λ₀)
        # all ℓ for fixed m
        @test @inferred(legendre(N, 0:ℓ, m, one(T))) isa Vector{T}
        @test @inferred(legendre!(N, Λ₁, LMAX, LMAX, one(T))) isa typeof(Λ₁)
        # all ℓ and m
        @test @inferred(legendre(N, 0:ℓ, 0:m, one(T))) isa Matrix{T}
        @test @inferred(legendre!(N, Λ₂, LMAX, LMAX, one(T))) isa typeof(Λ₂)
    end
end

@testset "Internal type promotion" begin
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

@testset "Broadcasting" begin
    x = 0.5
    z = range(-1, 1, length=10)
    Λ = zeros(length(z), LMAX + 1, LMAX + 1)

    # Single (l,m)
    @test λlm.(LMAX, LMAX, x) == λlm(LMAX, LMAX, x)
    @test λlm.(LMAX, LMAX, z) == λlm(LMAX, LMAX, z)
    # All l for fixed m
    @test λlm.(0:LMAX, LMAX, x) == λlm(0:LMAX, LMAX, x)
    @test λlm.(0:LMAX, LMAX, z) == λlm(0:LMAX, LMAX, z)
    # All l and m over multiple z
    @test λlm.(0:LMAX, 0:LMAX, x) == λlm(0:LMAX, 0:LMAX, x)
    @test λlm.(0:LMAX, 0:LMAX, z) == λlm(0:LMAX, 0:LMAX, z)

    # Test that the single (l,m) for scalar and vector arguments inferrs correctly, since
    # the return is branched on whether the argument is 0-dimensional or not.
    # Dot syntax is not a call expression, so have to test the broadcasted() call directly.
    @test @inferred(Base.broadcasted(legendre, λlm, LMAX, LMAX, x)) isa Float64
    @test @inferred(Base.broadcasted(legendre, λlm, LMAX, LMAX, z)) isa Vector{Float64}
end

@testset "Shape preservation ($(ndims(x))-dimensional array)" for x in (fill(1.0), ones(1), ones(1,1))
    @test ndims(λlm(LMAX, LMAX, x)) == ndims(x)
    @test ndims(λlm(0:LMAX, LMAX, x)) == ndims(x) + 1
    @test ndims(λlm(0:LMAX, 0:LMAX, x)) == ndims(x) + 2

    @test ndims(λlm.(LMAX, LMAX, x)) == ndims(x)
    @test ndims(λlm.(0:LMAX, LMAX, x)) == ndims(x) + 1
    @test ndims(λlm.(0:LMAX, 0:LMAX, x)) == ndims(x) + 2
end

@testset "Axes preservation (offset axes $axs)" for axs in ((-5:5,), (-5:5, -5:5))
    using OffsetArrays
    using Base: OneTo

    X = reshape(collect(range(-1, 1, length=prod(length.(axs)))), axs...)

    @test axes(λlm(LMAX, LMAX, X)) == axs
    @test axes(λlm(0:LMAX, LMAX, X)) == (axs..., OneTo(LMAX+1))
    @test axes(λlm(0:LMAX, 0:LMAX, X)) == (axs..., OneTo(LMAX+1), OneTo(LMAX+1))

    @test axes(λlm.(LMAX, LMAX, X)) == axs
    @test axes(λlm.(0:LMAX, LMAX, X)) == (axs..., OneTo(LMAX+1))
    @test axes(λlm.(0:LMAX, 0:LMAX, X)) == (axs..., OneTo(LMAX+1), OneTo(LMAX+1))

    Λ  = λlm(0:LMAX, 0:LMAX, parent(X)) # Normal axes
    Λ′ = λlm(0:LMAX, 0:LMAX, X)         # Offset axes
    @test parent(Λ′) == Λ # Equality of values
    @test_throws DimensionMismatch λlm!(Λ, LMAX, LMAX, X) # Mismatched axes
end

@testset "Preallocated work space" begin
    using AssociatedLegendrePolynomials: unsafe_legendre!
    x1 = 0.5
    xN = range(-1, 1, length=10)
    Λ1 = zeros(LMAX+1, LMAX+1)
    ΛN = zeros(length(xN), LMAX+1, LMAX+1)
    norm = LegendreUnitNorm()
    work1 = AssociatedLegendrePolynomials.Work(norm, Λ1, x1)
    workN = AssociatedLegendrePolynomials.Work(norm, ΛN, xN)
    # Check equality before allocations to ensure the methods have been compiled.
    @test unsafe_legendre!(norm, copy(Λ1), LMAX, LMAX, x1) ==
            unsafe_legendre!(work1, Λ1, LMAX, LMAX, x1)
    @test unsafe_legendre!(norm, copy(ΛN), LMAX, LMAX, xN) ==
            unsafe_legendre!(workN, ΛN, LMAX, LMAX, xN)
    # Normal version should allocate at least as much space as 8 working buffers
    @test 8sizeof(work1.z) ≤ @allocated unsafe_legendre!(norm, Λ1, LMAX, LMAX, x1)
    @test 8sizeof(workN.z) ≤ @allocated unsafe_legendre!(norm, ΛN, LMAX, LMAX, xN)
    # The pre-allocated buffer space should then permit zero allocation calls
    @test 0 == @allocated unsafe_legendre!(work1, Λ1, LMAX, LMAX, x1)
    @test 0 == @allocated unsafe_legendre!(workN, ΛN, LMAX, LMAX, xN)
end
