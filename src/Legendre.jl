"""
Collections of functions which compute the associated Legendre functions.

Based on implementation described in Limpanuparb and Milthorpe (2014)
*“Associated Legendre Polynomials and Spherical Harmonics Computation for
Chemistry Applications”* arXiv:1410.1748v1
"""
module Legendre

import Base: @boundscheck, @propagate_inbounds

include("scalar.jl")

# Public interfaces interface
export AbstractLegendreNorm
include("interface.jl")

# Specific normalizations
export LegendreUnitNorm,  LegendreSphereNorm,  LegendreNormCoeff,
       LegendreUnitCoeff, LegendreSphereCoeff
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
