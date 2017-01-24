module Harmonics
    using Base.Test

    include("harmonics/legendre.jl")

    function runtests()
        @testset "legendre" begin
            Legendre.runtests()
        end
    end
end

