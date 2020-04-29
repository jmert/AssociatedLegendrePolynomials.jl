using Test, TestSetExtensions
using Legendre # For doctests, include Legendre as a binding in Main

include("testsuite.jl")
using .TestSuite: LMAX, NumTypes
const Norms = (LegendreUnitNorm(), LegendreUnitCoeff{BigFloat}(LMAX))

function prettytime(t)
    v, u = t < 1e3 ? (t, "ns") :
           t < 1e6 ? (t/1e3, "μs") :
           t < 1e9 ? (t/1e6, "ms") :
                     (t/1e9, "s")
    return string(round(v, digits=3), " ", u)
end
macro include(file, desc)
    mod = gensym(first(splitext(file)))
    quote
        print($desc, ": ")
        t0 = time_ns()
        @testset $desc begin
            @eval module $mod
                using Test, Legendre
                import ..NumTypes
                include($file)
            end
        end
        printstyled("  ", prettytime(time_ns() - t0), "\n", color=:light_black)
    end
end

@testset ExtendedTestSet "Legendre" begin
    @include "norms.jl" "Normalizations"
    @include "coeffs.jl" "Precomputed coefficients"
    @include "implementation.jl" "Implementation details"
    @include "analytic.jl" "Analytic checks"
    @include "errors.jl" "Error checking"
    @include "scalar.jl" "Utility: scalar container"
    @include "doctests.jl" "Doctests"
end
