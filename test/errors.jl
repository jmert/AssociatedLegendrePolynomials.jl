@testset "Degree/order domain errors" begin
    LMAX = 10
    MMAX = 2
    Λ = Matrix{Float64}(undef, LMAX+1, MMAX+1)
    norms = (LegendreUnitNorm(), LegendreSphereNorm(),
             LegendreUnitCoeff{Float64}(LMAX, MMAX))

    @testset "ℓ < 0 ($(typeof(N))" for N in norms
        @test_throws DomainError legendre( N,    -1, 0, 0.5)
        @test_throws DomainError legendre!(N, Λ, -1, 0, 0.5)
    end
    @testset "ℓ > 0 & m < 0 ($(typeof(N))" for N in norms
        @test_throws DomainError legendre( N,    LMAX, -1, 0.5)
        @test_throws DomainError legendre!(N, Λ, LMAX, -1, 0.5)
    end
    @testset "m > ℓ ($(typeof(N))" for N in norms
        @test_throws DomainError legendre( N,    LMAX, LMAX+1, 0.5)
        @test_throws DomainError legendre!(N, Λ, LMAX, LMAX+1, 0.5)
    end
end

@testset "Degree/order range errors" begin
    LMAX = 5
    ctab = LegendreUnitCoeff{Float64}(LMAX)
    # Throws on invalid l ranges
    @test_throws ArgumentError Pl.(1:LMAX, 0.5)
    @test_throws ArgumentError Plm.(1:LMAX, 0, 0.5)
    @test_throws ArgumentError λlm.(1:LMAX, 0, 0.5)
    @test_throws ArgumentError ctab.(1:LMAX, 0, 0.5)
    @test_throws ArgumentError legendre.(ctab, 1:LMAX, 0, 0.5)
    # Throws on invalid m ranges
    @test_throws ArgumentError Plm.(0:LMAX, 1:LMAX, 0.5)
    @test_throws ArgumentError λlm.(0:LMAX, 1:LMAX, 0.5)
    @test_throws ArgumentError ctab.(0:LMAX, 1:LMAX, 0.5)
    @test_throws ArgumentError legendre.(ctab, 0:LMAX, 1:LMAX, 0.5)
end

@testset "Output array bounds checking" begin
    LMAX = 2

    # Scalar argument
    λ  = Vector{Float64}(undef, LMAX)
    Λ₁ = Matrix{Float64}(undef, LMAX, LMAX+1)
    Λ₂ = Matrix{Float64}(undef, LMAX+1, LMAX)
    # Insufficient dim 1 (ℓ)
    @test_throws DimensionMismatch Plm!(λ, LMAX, 0, 0.5)
    @test_throws DimensionMismatch Plm!(Λ₁, LMAX, LMAX, 0.5)
    # Insufficient dim 2 (m)
    @test_throws DimensionMismatch Plm!(Λ₂, LMAX, LMAX, 0.5)

    # Vector argument
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

@testset "Precompute coefficient lmax/mmax" begin
    LMAX, MMAX = 2, 0
    ctab = LegendreUnitCoeff{Float64}(LMAX, MMAX)
    @test_throws BoundsError ctab(LMAX+1, 0, 0.5)
    @test_throws BoundsError ctab(LMAX, MMAX+1, 0.5)
end

