```@meta
CurrentModule = Legendre
```
# API Reference

## Functions

```@docs
legendre
legendre!
Nlm
```

## Normalizations

```@docs
AbstractLegendreNorm
LegendreUnitNorm
LegendreFourPiNorm
LegendreOrthoNorm
LegendreSphereNorm
LegendreNormCoeff
```

## Aliases

The following functors are constant aliases to the underlying normalization types which have
been made callable via an appropriate call overload.

```@docs
Plm
λlm
Plm!
λlm!
```

There are also aliases for pre-computed coefficients of the provided normalizations.
```@docs
LegendreUnitCoeff
LegendreFourPiCoeff
LegendreOrthoCoeff
LegendreSphereCoeff
```

## Normalization Interface

The following functions are unexported but considered a public API for interacting with
and defining additional normalizations.

```@docs
initcond
coeff_μ
coeff_ν
coeff_α
coeff_β
boundscheck_hook
```
