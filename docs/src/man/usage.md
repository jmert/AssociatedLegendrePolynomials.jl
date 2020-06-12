# [Usage](@id usage)
```@meta
DocTestFilters = Regex[
        r"Ptr{0x[0-9a-f]+}",
        r"[0-9\.]+ seconds( \(.*\))?",
        ]
```

```@contents
Pages = ["usage.md"]
Depth = 2
```

## Calculating scalar values

At its simplest, the associated Legendre polynomial ``P_ℓ^m(x)`` is computed by calling
[`Legendre.Plm`](@ref). For example, to compute ``P_2^1(0.5)``,
```jldoctest PlmUsage
julia> using Legendre

julia> Plm(2, 1, 0.5)
-1.299038105676658
```
In the context of CMB analysis, a common use of the associated Legendre polynomials is to
compute the spherical harmonics ``Y_{ℓm}(θ,ϕ)``:
```math
\begin{align}
    \begin{aligned}
    Y_{ℓm}(θ,ϕ) &≡ N_ℓ^m P_ℓ^m(\cos θ) e^{imϕ} \\
    &\text{where } N_ℓ^m ≡ \sqrt{\frac{2ℓ+1}{4π} \frac{(ℓ-m)!}{(ℓ+m)!}}
    \end{aligned}
\end{align}
```
The function [`Legendre.Nlm`](@ref) calculates the normalization factor ``N_ℓ^m``:
```jldoctest PlmUsage
julia> Nlm(2, 0)
0.6307831305050401

julia> Nlm(2, 0) * Plm(2, 0, 0.5)
-0.07884789131313001
```

An important fact about the associated Legendre polynomials is that for
``m > 0``, ``P_ℓ^m(x)`` diverges to ``∞`` as ``ℓ → ∞`` [^1].
For even moderately large pairs of ``(ℓ,m)``, numerical underflow and overflow make
computing the spherical harmonics impossible this way:
```jldoctest PlmUsage
julia> n = Nlm(157, 150)      # Underflows
0.0

julia> p = Plm(157, 150, 0.5) # Overflows
Inf

julia> n * p                  # Undefined
NaN
```

One way around this would be to just use extended precision arithmetic
```jldoctest PlmUsage
julia> n = Nlm(BigFloat, 157, 150)
4.14800666209481424285411223457923933542541063872695815968861285171699012214351e-314

julia> p = Plm(157, 150, big"0.5")
4.768286486602206390406601862422168575170463348990958242752608686436785229641823e+308

julia> Float64(n * p)
1.9778884113202627e-5
```
but at the expense of much more computationally expensive calculations.

An alternative way forward is to directly calculate the spherical harmonic normalized
associated Legendre polynomials ``λ_ℓ^m(x)`` so that the spherical harmonics are
defined as
```math
\begin{align}
    \begin{aligned}
    Y_{ℓm}(θ,ϕ) &= λ_ℓ^m(\cos θ) e^{imϕ} \\
    & \text{where } λ_ℓ^m(x) ≡ N_ℓ^m P_ℓ^m(x)
    \end{aligned}
\end{align}
```
[`Legendre.λlm`](@ref) implements this scheme and avoids the under/overflow of
computing the normalization separately from the function:
```jldoctest PlmUsage
julia> λlm(157, 150, 0.5)
1.977888411320263e-5
```

!!! note
    We are not just limited to efficient and numerically stable computation of
    ``λ_ℓ^m(x)``; the package supports arbitrary normalizations.  For further
    information on implementing custom Legendre normalizations, see the [Custom
    normalizations](@ref customnorm) section.

## Calculating multiple degrees/orders

Because calculating a particular Legendre polynomial value is the end result of running
a recurrence relation, looping evaluation of ``P_ℓ^m(x)`` for all ``ℓ`` is inefficient and
redoes a lot of work:
```jldoctest PlmUsage
julia> @time [l < 2 ? 0.0 : λlm(l, 2, 0.5) for l in 2:700];
  0.039210 seconds (71.21 k allocations: 3.285 MiB)
```
It's far more efficient to accumulate the intermediate terms while running the recurrence
relations.
Using a `UnitRange` as the input degree causes the functions to allocate and fill the
vector with all polynomials values:
```jldoctest PlmUsage
julia> λ = @time λlm(0:700, 2, 0.5);
  0.000012 seconds (6 allocations: 5.703 KiB)
```
On my machine, this is roughly 3000 times faster!

Likewise, calculating the [lower triangular] matrix of values for some ``x`` over all
degrees ``ℓ ∈ [0,ℓ_\mathrm{max}]`` and all orders ``m ∈ [0,ℓ]`` is done by also specifying
the orders as a `UnitRange`.
```jldoctest PlmUsage
julia> Λ = @time λlm(0:700, 0:700, 0.5);
  0.002980 seconds (7 allocations: 3.749 MiB)

julia> Λ[:,3] == λ   # N.B. 1-based indexing of the array!
true
```

!!! note
    There are two things in particular to remember with the range-based calls:
    1. The ranges must start at 0, otherwise an `ArgumentError` will be thrown.
    2. Calculating a range of orders ``m`` for a fixed degree ``\ell`` is not supported;
       to calculate multiple orders requires the output matrix be at least square (but
       may "tall and skinny" with ``\ell_\mathrm{max} > m_\mathrm{max}`` if desired).

It is also more efficient to operate upon an array of arguments ``x`` than to loop over
them one-by-one, so the functions also accept the input argument `x` as an array of
any shape.

For a specific degree and order, the output array will have the same shape as the
argument:
```jldoctest PlmUsage
julia> λlm(2, 0, reshape(range(0, 1, length=4), 2, 2))
2×2 Array{Float64,2}:
 -0.315392  0.105131
 -0.210261  0.630783
```
Then adding a range of degrees increases the dimensionality by 1, with the trailing
dimension being over ``\ell``,
```jldoctest PlmUsage
julia> λlm(0:2, 0, reshape(range(0, 1, length=4), 2, 2))
2×2×3 Array{Float64,3}:
[:, :, 1] =
 0.282095  0.282095
 0.282095  0.282095

[:, :, 2] =
 0.0       0.325735
 0.162868  0.488603

[:, :, 3] =
 -0.315392  0.105131
 -0.210261  0.630783
```
and a further extra dimension is added for a range over orders ``m``.

## In-place calculations

Both of `Plm` and `λlm` also have in-place modifying counterparts,
[`Plm!`](@ref Plm!(::Any, ::Integer, ::Integer, ::Any)) and
[`λlm!`](@ref λlm!(::Any, ::Integer, ::Integer, ::Any)) respectively,
which fill an appropriately sized vector for a specified ``ℓ_\mathrm{max}`` and
``m_\mathrm{max}``.
Instead of using integer or range arguments, whether to calculate a value for a
single degree/order, a range of degrees for fixed order, or for all degrees and orders
is inferred based on the dimensionality of the output array.

For example, to calculate the single value ``\lambda_{700}^{200}(0.5)``, provide a
0-dimensional output array (to match the 0-dimensionality of the scalar `0.5`)
```jldoctest PlmUsage
julia> λlm!(fill(NaN), 700, 2, 0.5)
0-dimensional Array{Float64,0}:
0.24148976866924293
```
and filling a vector or matrix instead calculates all degrees up to the given maximum
degree/order as appropriate:
```jldoctest PlmUsage
julia> λlm!(λ, 700, 2, 0.5) == λlm(0:700, 2, 0.5)
true

julia> λlm!(Λ, 700, 700, 0.5) == λlm(0:700, 0:700, 0.5)
true
```

The in-place interface accepts input arguments `x` of any shape as well, with the output
array `Λ` having to have between 0 and 2 more dimensions than `x`, where the leading
dimensions of the input and output arrays have the same axes, and the trailing dimensions
are sized appropriate for the number of degrees/orders to be calculated.

## Precomputed recursion factors

A final trick to accelerating calculation of any normalization of the associated
Legendre polynomials is to pre-compute the appropriate recursion relation coefficients.

At a low level, `Plm`/`Plm!` and `λlm`/`λlm!` are simple wrappers around the general
[`legendre`](@ref)/[`legendre!`](@ref) functions.
The trait type [`LegendreUnitNorm`](@ref) dispatches internal functions to compute
``P_ℓ^m(x)``:
```jldoctest PlmUsage
julia> legendre(LegendreUnitNorm(), 5, 2, 0.5) == Plm(5, 2, 0.5)
true
```
 and [`LegendreSphereNorm`](@ref) does the same for ``λ_ℓ^m(x)``:
```jldoctest PlmUsage
julia> legendre(LegendreSphereNorm(), 5, 2, 0.5) == λlm(5, 2, 0.5)
true
```

The type [`LegendreNormCoeff`](@ref) stores the coefficients for a particular
normalization (and value type) so that the coefficients must only be calculated once.
Aliases for the unit and spherical normalizations are provided by default,
[`LegendreUnitCoeff`](@ref) and [`LegendreSphereCoeff`](@ref) respectively.
```jldoctest PlmUsage
julia> coeff = LegendreSphereCoeff{Float64}(700);

julia> legendre(coeff, 5, 2, 0.5)
-0.15888479843070935
```

!!! warning "Performance Note"

    Choosing whether to use the pre-computed coefficients or not should be guided by
    benchmarking and performance profiling.
    Modern processors can perform many floating point operations in the time it takes
    to load the coefficients from memory, so depending on the complexity of the
    normalization, you may actually achieve better performance by recomputing the
    recursion coefficients on demand.

Notice that due to its flexibility, `legendre!` requires explicit `lmax` and `mmax`
arguments even though the `LegendreNormCoeff` has a `lmax` and `mmax` set during
construction.
This allows us to pass both a coefficient cache and output array which are larger than the
computed set of coefficients.
For example, the output matrix and cache used above each support computing the Legendre
polynomials up to ``\ell = 700``, but if we only need ``\ell \le 2``, we can avoid
computing terms beyond our required problem size.
```jldoctest PlmUsage
julia> fill!(Λ, 0);

julia> legendre!(coeff, Λ, 2, 2, 0.5);

julia> Λ[1:5, 1:5]
5×5 Array{Float64,2}:
  0.282095    0.0       0.0       0.0  0.0
  0.244301   -0.299207  0.0       0.0  0.0
 -0.0788479  -0.334523  0.289706  0.0  0.0
  0.0         0.0       0.0       0.0  0.0
  0.0         0.0       0.0       0.0  0.0
```

---

## Footnotes

[^1]:
    Specifically, the envelope of ``P_ℓ^m(x)`` which bounds the local extrema
    for all values of ``x`` can be shown to be
    ```math
        \left| P_ℓ^m(\cos θ) \right| ≤ \frac{Γ(ℓ+m+1)}{Γ(ℓ+\frac{3}{2})}
            \left( \frac{2}{π \sin θ} \right)^{1/2}
    ```
    (see Eq. 8.10.7 (p336) of Abramowitz and Stegun, “Handbook of Mathematical
    Functions” 10th printing (1972)).
    For fixed ``m`` and any ``x``, we take the asymptotic limit as ``ℓ → ∞`` and simplify
    ``Γ(z)`` via Stirling's approximation to get the scaling of the associated Legendre
    polynomial envelope
    ```math
        \DeclareMathOperator*{\env}{env}
        \env_{ℓ→∞}\left( P_ℓ^m \right) ∝ ℓ^{m - 1/2} \text{ .}
    ```
    In contrast, the normalization factor ``N_ℓ^m`` scales as ``ℓ^{1/2 - m}``,
    exactly canceling the scaling of ``\env\left(P_ℓ^m\right)``, so overall the spherical
    harmonic normalized Legendre polynomials ``λ_ℓ^m(x)`` asymptote to some constant
    envelope:
    ```math
        \env_{ℓ→∞} \left( λ_ℓ^m \right) ∝ ℓ^0 = \text{constant .}
    ```
