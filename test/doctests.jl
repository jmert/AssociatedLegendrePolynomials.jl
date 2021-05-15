using Documenter, Logging

# Outputs are Julia version specific
if VERSION >= v"1.6"
    # Disable Documeter's Info logging
    oldlvl = Logging.min_enabled_level(current_logger())
    disable_logging(Logging.Info)
    try
        DocMeta.setdocmeta!(Legendre, :DocTestSetup, :(using Legendre); recursive=true)
        doctest(Legendre, testset="Doc Tests")
    finally
        disable_logging(oldlvl - 1)
    end
end
