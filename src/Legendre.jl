module Legendre

import Base: @boundscheck, @propagate_inbounds

include("scalar.jl")

# Public interfaces interface
export AbstractLegendreNorm
include("interface.jl")

# Specific normalizations
export LegendreUnitNorm,  LegendreOrthoNorm,  LegendreSphereNorm,  LegendreNormCoeff,
       LegendreUnitCoeff, LegendreOrthoCoeff, LegendreSphereCoeff
include("norm_unit.jl")
include("norm_sphere.jl")
include("norm_table.jl")

export legendre, legendre!
include("calculation.jl")

# Other functionality
export Plm, Plm!, Nlm, λlm, λlm!
include("aliases.jl")
include("broadcasting.jl")

include("precompile.jl")
_precompile_()
end # module Legendre
