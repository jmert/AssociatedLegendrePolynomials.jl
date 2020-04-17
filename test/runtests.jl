using Test, Logging
using Legendre
const NumTypes = (Float32, Float64, BigFloat)

const TESTLIST = [
    "legendre" => "Legendre",
   ]

@testset "Legendre" begin
    @testset "$desc" for (id,desc) in TESTLIST
        modpath = joinpath(dirname(@__FILE__), "$(id).jl")
        # Include the file and have it processed within this module
        print("running $desc tests... ")
        t0 = time_ns()
        include(modpath)
        t1 = time_ns()
        println( (t1-t0)/1e9, " seconds")
    end

    # Disable Documeter's Info logging
    oldlvl = Logging.min_enabled_level(current_logger())
    disable_logging(Logging.Info)
    try
        using Documenter
        DocMeta.setdocmeta!(Legendre, :DocTestSetup, :(using Legendre); recursive=true)
        print("running Doc tests... ")
        t0 = time_ns()
        doctest(Legendre, testset="Doc Tests")
        t1 = time_ns()
        println( (t1-t0)/1e9, " seconds")
    finally
        disable_logging(oldlvl - 1)
    end
end
