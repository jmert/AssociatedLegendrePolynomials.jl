import ..TestSuite

@testset "Unit normalization" begin
    TestSuite.runtests(LegendreUnitNorm())
end
@testset "Spherical normalization" begin
    TestSuite.runtests(LegendreSphereNorm())
end
