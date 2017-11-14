# Legendre Polynomials

```@contents
Pages = ["legendre.md"]
Depth = 2
```

The [`Legendre`](@ref) module implementation has been largely based on the approach of
[Limpanuparb and Milthorpe (2014)](@ref bib-legendre).

## [Definition and Properties](@id legendre_defn)

The associated Legendre polynomials ``P_ℓ^m(x)`` are the solution to the differential
equation
```math
\begin{align}
    (1-x^2) \frac{d^2}{dx^2}P_ℓ^m(x) - 2x \frac{d}{dx}P_ℓ^m(x) + \left[ ℓ(ℓ+1) -
        \frac{m^2}{1-x^2} \right] P_ℓ^m(x) = 0
\end{align}
```
which arises as the colatitude ``θ`` part of solving Laplace's equation
``∇^2 ψ + λψ = 0`` in spherical coordinates (where ``x = \cos(θ)``).

There are several different conventions used to define ``P_ℓ^m`` that provide
different properties, but the convention used here is typical of quantum
mechanics and obeys the following properties:

* Solutions only exist for integer ``ℓ`` and ``m``, where ``ℓ ≤ 0`` and ``|m| ≤ ℓ``.

* The associated Legendre functions are normalized such that ``P_0^0`` is unity and have
  orthogonality conditions,
  ```math
  \begin{align}
      \int_{-1}^1 P_ℓ^m(x) P_{ℓ'}^{m}(x)\,\mathrm{d}x
          = \frac{2}{2ℓ+1} \frac{(ℓ+m)!}{(ℓ-m)!}
          \delta_{ℓℓ'}
  \end{align}
  ```
  for constant ``m`` and
  ```math
  \begin{align}
      \int_{-1}^1 \frac{P_ℓ^m(x) P_{ℓ}^{m'}(x)}{1-x^2}\,\mathrm{d}x
          = \frac{1}{m} \frac{(ℓ+m)!}{(ℓ-m)!} \delta_{mm'}
  \end{align}
  ```
  for constant ``ℓ``, where ``δ`` is the Kronecker delta.

* The phase convention for the Legendre functions is chosen such that the negative orders
  are related to positive orders according to,
  ```math
  \begin{align}
      P_ℓ^{-m}(x) = (-1)^m \frac{(ℓ-m)!}{(ℓ+m)!} P_ℓ^m(x)
  \end{align}
  ```

* The Legendre functions can be enumerated for non-negative ``m`` using the three
  following recursion relations (given the initial condition ``P_0^0(x)``):
  ```math
  \begin{align}
      (ℓ - m + 1)P_{ℓ+1}^m(x) &= (2ℓ+1)xP_ℓ^m(x) - (ℓ+m)P_{ℓ-1}^m(x)
      \label{eqn:std_rr_2term}
      \\
      P_{ℓ+1}^{ℓ+1}(x) &= -(2ℓ+1)\sqrt{1-x^2} P_ℓ^ℓ(x)
      \label{eqn:std_rr_1term_lm}
      \\
      P_{ℓ+1}^ℓ(x) &= x(2ℓ+1)P_ℓ^ℓ(x)
      \label{eqn:std_rr_1term_l}
  \end{align}
  ```

## [Usage](@id legendre_usage)

### Calculating scalar values

At the most basic, the associated Legendre polynomial ``P_ℓ^m(x)`` is computed by calling
[`CMB.Legendre.Plm`](@ref). For example, to compute ``P_2^1(0.5)``,
```jldoctest PlmUsage
julia> using CMB.Legendre

julia> Plm(2, 1, 0.5)
-1.299038105676658
```
When ``m = 0`` and only the Legendre polynomial ``P_ℓ(x)`` is needed,
[`CMB.Legendre.Pl`](@ref) can be used instead:
```jldoctest PlmUsage
julia> Plm(2, 0, 0.5)
-0.125

julia> Pl(2, 0.5)
-0.125

julia> Pl(2, 0.5) == Plm(2, 0, 0.5)
true
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
The function [`CMB.Legendre.Nlm`](@ref) calculates the normalization factor ``N_ℓ^m``:
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
4.768286486602206390406601862422168575170463348990958242752608686436785229641202e+308

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
[`CMB.Legendre.λlm`](@ref) implements this scheme and avoids the under/overflow of
computing the normalization separately from the function:
```jldoctest PlmUsage
julia> λlm(157, 150, 0.5)
1.977888411320241e-5
```

!!! note
    We are not just limited to efficient and numerically stable computation of
    ``λ_ℓ^m(x)``; the package supports arbitrary normalizations.  For further
    information on implementing custom Legendre normalizations, see the [Custom
    normalizations](@ref legendre_customnorm) section.

### Calculating all values up to a given ``ℓ_\mathrm{max}``

Because calculating a particular Legendre polynomial value is the end result of running
a recurrence relation, using Julia's dot broadcasting to compute ``P_ℓ^m(x)``
for all ``ℓ`` is inefficient and redoes a lot of work:
```julia
julia> Λ = zeros(701);

julia> @time Λ[3:701] .= λlm.(2:700, 2, 0.5);
  0.042107 seconds (4.61 k allocations: 257.940 KiB)
```
It's far more efficient to incrementally calculate the ``ℓ+1`` term directly from the
``ℓ`` term.
Both of `Plm` and `λlm` have modifying counterparts,
[`Plm!`](@ref Plm!(::AbstractVector, ::Integer, ::Integer, ::Real)) and
[`λlm!`](@ref λlm!(::AbstractVector, ::Integer, ::Integer, ::Real)) respectively,
which fill an appropriately sized vector for a specified ``ℓ_\mathrm{max}``.
```julia
julia> @time λlm!(Λ, 700, 2, 0.5);
  0.000036 seconds (4 allocations: 160 bytes)
```
On my machine, this ends up being roughly 1000 times faster!

If all Legendre polynomial values for some ``x`` over all
``ℓ ∈ [0,ℓ_\mathrm{max}]`` and ``m ∈ [0,ℓ]`` are required, there are also methods of
[`Plm!`](@ref Plm!(::AbstractMatrix, ::Integer, ::Real)) and
[`λlm!`](@ref λlm!(::AbstractMatrix, ::Integer, ::Real))
which fill the entire [lower triangular] matrix of values:
```jldoctest PlmUsage
julia> Λ = zeros(701, 701);

julia> λlm!(Λ, 700, 0.5);

julia> Λ[701,3] == λlm(700, 2, 0.5)   # N.B. 1-based indexing of the array!
true
```

### Precomputed recursion factors

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
```jldoctest PlmUsage
julia> coeff = LegendreNormCoeff{LegendreSphereNorm, Float64}(700);

julia> legendre(coeff, 5, 2, 0.5)
-0.15888479843070935
```
On my machine, this results in a further ~50% decrease in computation time compared to
`λlm!`:
```julia
julia> @time legendre!(coeff, Λ, 700, 2, 0.5);
  0.000020 seconds (4 allocations: 160 bytes)
```

## [Custom normalizations](@id legendre_customnorm)

---

### Footnotes

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
