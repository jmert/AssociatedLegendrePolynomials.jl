import ..TestSuite

@testset "Precomputed coefficients" begin
    TestSuite.runtests(LegendreUnitCoeff{Float64}(TestSuite.LMAX))
    TestSuite.runtests(LegendreUnitCoeff{Float32}(TestSuite.LMAX))
    TestSuite.runtests(LegendreUnitCoeff{BigFloat}(TestSuite.LMAX))

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
end

