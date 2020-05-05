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
    Y_{ℓm}(θ,ϕ) &≡ (-1)^m N_ℓ^m P_ℓ^m(\cos θ) e^{imϕ} \\
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
    Y_{ℓm}(θ,ϕ) &= (-1)^m λ_ℓ^m(\cos θ) e^{imϕ} \\
    & \text{where } λ_ℓ^m(x) ≡ N_ℓ^m P_ℓ^m(x)
    \end{aligned}
\end{align}
```
[`Legendre.λlm`](@ref) implements this scheme and avoids the under/overflow of
computing the normalization separately from the function:
```jldoctest PlmUsage
julia> λlm(157, 150, 0.5)
1.977888411320258e-5
```

!!! note
    We are not just limited to efficient and numerically stable computation of
    ``λ_ℓ^m(x)``; the package supports arbitrary normalizations.  For further
    information on implementing custom Legendre normalizations, see the [Custom
    normalizations](@ref customnorm) section.

## Calculating all values up to a given ``ℓ_\mathrm{max}``

Because calculating a particular Legendre polynomial value is the end result of running
a recurrence relation, looping evaluation of ``P_ℓ^m(x)`` for all ``ℓ`` is inefficient and
redoes a lot of work:
```jldoctest PlmUsage
julia> λ = zeros(701);

julia> @time λ[3:701] = [λlm(l, 2, 0.5) for l in 2:700];
  0.063346 seconds (56.42 k allocations: 2.539 MiB)
```
It's far more efficient to accumulate the intermediate terms while running the recurrence
relations.
Both of `Plm` and `λlm` have modifying counterparts,
[`Plm!`](@ref Plm!(::Any, ::Integer, ::Integer, ::Any)) and
[`λlm!`](@ref λlm!(::Any, ::Integer, ::Integer, ::Any)) respectively,
which fill an appropriately sized vector for a specified ``ℓ_\mathrm{max}``.
```jldoctest PlmUsage
julia> @time λlm!(λ, 700, 2, 0.5);
  0.000162 seconds (14 allocations: 320 bytes)
```
On my machine, this ends up being roughly 400 times faster!

If all Legendre polynomial values for some ``x`` over all
``ℓ ∈ [0,ℓ_\mathrm{max}]`` and ``m ∈ [0,ℓ]`` are required, instead supply an output
matrix into which the lower triangle of values is filled:
```jldoctest PlmUsage
julia> Λ = zeros(701, 701);

julia> λlm!(Λ, 700, 700, 0.5);

julia> Λ[701,3] == λlm(700, 2, 0.5)   # N.B. 1-based indexing of the array!
true
```

## Broadcasting over multiple arguments

The Legendre polynomials can be evaluated over multiple arguments ``x`` as well by
using Julia's standard broadcasting syntax:
```jldoctest PlmUsage
julia> λlm.(2, 0, range(-1.0, 1.0, length=5))
5-element Array{Float64,1}:
  0.63078313050504
 -0.07884789131313
 -0.31539156525252
 -0.07884789131313
  0.63078313050504
```
Broadcasting has been specialized for calls to `Pl`, `Plm`, and `λlm` to avoid the
overhead inherent in calling the scalar functions multiple times:
```jldoctest PlmUsage
julia> z = range(-1.0, 1.0, length=500);

julia> @time [λlm(2, 0, i) for i in z];
  0.054037 seconds (54.23 k allocations: 2.447 MiB)

julia> @time λlm.(2, 0, z);
  0.000032 seconds (13 allocations: 36.719 KiB)
```
In fact, the shape of `z` is preserved, so any matrix shape can be used:
```jldoctest PlmUsage; setup=(using Random; Random.seed!(0))
julia> λlm.(2, 0, rand(3,3,3))
3×3×3 Array{Float64,3}:
[:, :, 1] =
  0.326489  -0.285639  -0.313698
  0.46875   -0.241804  -0.310982
 -0.289767  -0.276217  -0.191519

[:, :, 2] =
  0.580778   -0.251413   0.091097
  0.0093121   0.468216  -0.00159624
 -0.0402128  -0.288992   0.397937

[:, :, 3] =
  0.57083   -0.311711  -0.313631
  0.242235  -0.197404  -0.247441
 -0.107      0.242107  -0.311164
```

Obtaining the Legendre polynomials over multiple ``\ell`` and/or ``m`` values for many
arguments can be done via broadcasting as well.
The degree `l` must be a `UnitRange` starting at zero, and `m` may be either a
scalar integer (to calculate all ``\ell`` for a fixed ``m``) or a `UnitRange`
starting at zero as well.
For example, to compute the ``P_\ell^0(z)`` and ``P_\ell^2(z)`` coefficients to an
``\ell_{\mathrm{max}} = 700``, one could execute
```jldoctest PlmUsage
julia> summary(Plm.(0:700, 0:2, z))
"500×701×3 Array{Float64,3}"
```
The output array will have between 0 and 2 more dimensions more than the dimensionality
of the input arguments depending on the calling convention.
For scalar values of `l` and `m`, the output will be the same shape as `z` with no
extra trailing dimensions.
If instead `l` is a `UnitRange`, the output dimensionality increases by one, and the
trailing dimension runs over the degrees ``\ell``;
switching to `m` a `UnitRange` as well, the output dimensionality is two greater than
`z`, with the penultimate and final dimensions running over ``\ell`` and ``m``,
respectively.

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
