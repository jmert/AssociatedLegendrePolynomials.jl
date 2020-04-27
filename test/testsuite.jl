module TestSuite
    using Test, Legendre
    export runtests

    function runtests(norm)
        interface(norm)
    end

    function interface(norm)
        @testset "Legendre Normalization ($norm)" begin
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
end
