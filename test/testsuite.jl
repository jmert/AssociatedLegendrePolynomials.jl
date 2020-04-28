module TestSuite
    using Test, Legendre
    const LMAX = 5
    const NumTypes = (Float32, Float64, BigFloat)

    function runtests(norm)
        @testset "Normalization $(typeof(norm))" begin
            interface(norm)
            promotion(norm)
        end
    end

    function interface(norm)
        # Test all the aspects of the interface that each normalization type must
        # implement
        @testset "Legendre Interface" begin
            @test norm isa AbstractLegendreNorm
            # Initial condition
            @test isfinite(Legendre.Plm_00(norm, Float64))
            # 1-term recurrence coefficients
            @test isfinite(Legendre.Plm_μ(norm, Float64, 1))
            @test isfinite(Legendre.Plm_ν(norm, Float64, 1))
            # 2-term recurrence coefficients
            @test isfinite(Legendre.Plm_α(norm, Float64, 1, 0))
            @test isfinite(Legendre.Plm_β(norm, Float64, 1, 0))
        end
    end

    function promotion(norm)
        Λ₀ = fill(0.0)
        Λ₁ = zeros(LMAX+1)
        Λ₂ = zeros(LMAX+1, LMAX+1)
        @testset "Argument eltype $T" for T in NumTypes
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
                @test @inferred(legendre(norm, ℓ, m, one(T))) isa T
                # Also must handle single value, all ℓ for fixed m, and all (ℓ,m) up to
                # (LMAX,LMAX)
                @test @inferred(legendre!(norm, Λ₀, LMAX, LMAX, one(T))) isa typeof(Λ₀)
                @test @inferred(legendre!(norm, Λ₁, LMAX, LMAX, one(T))) isa typeof(Λ₁)
                @test @inferred(legendre!(norm, Λ₂, LMAX, LMAX, one(T))) isa typeof(Λ₂)
            end
        end
    end
end
