using Test, TestSetExtensions, Logging
const NumTypes = (Float32, Float64, BigFloat)

macro include(file, desc)
    mod = Symbol(splitext(file))
    quote
        print($desc, ": ")
        @eval module $mod
            using Test, Legendre
            import ..NumTypes
            @testset $desc begin
                Base.include($mod, $file)
            end
        end
        println()
    end
end

@testset ExtendedTestSet "Legendre" begin
    @include "scalar.jl" "Broadcastable scalar"
    @include "legendre.jl" "Legendre"

    print("Doc tests: ")
    # Disable Documeter's Info logging
    oldlvl = Logging.min_enabled_level(current_logger())
    disable_logging(Logging.Info)
    try
        using Documenter, Legendre
        DocMeta.setdocmeta!(Legendre, :DocTestSetup, :(using Legendre); recursive=true)
        doctest(Legendre, testset="Doc Tests")
    finally
        disable_logging(oldlvl - 1)
    end
end
