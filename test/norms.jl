import ..TestSuite

@testset "Unit normalization" begin
    TestSuite.runtests(LegendreUnitNorm())
end
@testset "Orthonormal normalization" begin
    TestSuite.runtests(LegendreOrthoNorm())
end
@testset "Spherical normalization" begin
    TestSuite.runtests(LegendreSphereNorm())
end
