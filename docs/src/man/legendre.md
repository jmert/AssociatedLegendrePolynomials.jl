# Legendre Polynomials
```@meta
DocTestFilters = Regex[
        r"Ptr{0x[0-9a-f]+}",
        r"[0-9\.]+ seconds( \(.*\))?",
        ]
```

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
[`CMB.Legendre.λlm`](@ref) implements this scheme and avoids the under/overflow of
computing the normalization separately from the function:
```jldoctest PlmUsage
julia> λlm(157, 150, 0.5)
1.977888411320258e-5
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
```jldoctest PlmUsage
julia> λ = zeros(701);

julia> @time λ[3:701] .= λlm.(2:700, 2, 0.5);
  0.042107 seconds (4.61 k allocations: 257.940 KiB)
```
It's far more efficient to incrementally calculate the ``ℓ+1`` term directly from the
``ℓ`` term.
Both of `Plm` and `λlm` have modifying counterparts,
[`Plm!`](@ref Plm!(::AbstractVector, ::Integer, ::Integer, ::Real)) and
[`λlm!`](@ref λlm!(::AbstractVector, ::Integer, ::Integer, ::Real)) respectively,
which fill an appropriately sized vector for a specified ``ℓ_\mathrm{max}``.
```jldoctest PlmUsage
julia> @time λlm!(λ, 700, 2, 0.5);
  0.000036 seconds (4 allocations: 160 bytes)
```
On my machine, this ends up being roughly 1000 times faster!

If all Legendre polynomial values for some ``x`` over all
``ℓ ∈ [0,ℓ_\mathrm{max}]`` and ``m ∈ [0,ℓ]`` are required, there are also methods of
[`Plm!`](@ref Plm!(::AbstractMatrix, ::Integer, ::Integer, ::Real)) and
[`λlm!`](@ref λlm!(::AbstractMatrix, ::Integer, ::Integer, ::Real))
which fill the entire [lower triangular] matrix of values:
```jldoctest PlmUsage
julia> Λ = zeros(701, 701);

julia> λlm!(Λ, 700, 700, 0.5);

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
Aliases for the unit and spherical normalizations are provided by default,
[`LegendreUnitCoeff`](@ref) and [`LegendreSphereCoeff`](@ref) respectively.
```jldoctest PlmUsage
julia> coeff = LegendreSphereCoeff{Float64}(700);

julia> legendre(coeff, 5, 2, 0.5)
-0.15888479843070935
```
On my machine, this results in a further ~50% decrease in computation time compared to
`λlm!`:
```jldoctest PlmUsage
julia> @time legendre!(coeff, λ, 700, 2, 0.5);
  0.000020 seconds
```

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
In most situations, though, it'll probably be most convenient to use the functor interface
attached to the coefficient cache object which assumes the `lmax` it was constructed
with.
The coefficient table itself is callable with forms similar to `legendre` and `legendre!`
except that the `norm` and `lmax` arguments are implicit/not necessary.
```jldoctest PlmUsage
julia> coeff(20, 0.5)    # == legendre(coeff, 20, 0.5)
-0.08734916334699527

julia> coeff(20, 2, 0.5) # == legendre(coeff, 20, 2, 0.5)
0.10617507806374693

julia> leg! = coeff;    # alias to clarify that leg! modifies

julia> leg!(λ, 2, 0.5); # same as legendre!(coeff, λ, size(coeff.α, 1) - 1, 2, 0.5)

julia> leg!(Λ, 0.5);    # same as legendre!(coeff, Λ, (size(coeff.α) .- 1)..., 0.5)
```

## [Custom normalizations](@id legendre_customnorm)
`CMB.Legendre` provides the standard and spherical harmonic normalizations by default, but
arbitrary normalizations are also supported.
The mile-high overview is that the initial condition and recurrence relation (r.r.)
coefficients are all methods which dispatch on a normalization trait type, so a new
normalization is added by simply extending appropriate types and methods.
The following table lists all of the types to extend and method specialization to
implement.

| Interfaces to extend/implement              | Brief description                                                              |
|:------------------------------------------- |:------------------------------------------------------------------------------ |
| [`CMB.Legendre.AbstractLegendreNorm`](@ref) | Supertype of normalization trait types                                         |
| [`CMB.Legendre.Plm_00()`](@ref)             | Value of ``N_0^0 P_0^0(x)`` for the given normalization                        |
| [`CMB.Legendre.Plm_μ()`](@ref)              | Coefficient ``μ_m`` for the 1-term r.r. boosting ``ℓ → ℓ+1`` and ``m → m+1``   |
| [`CMB.Legendre.Plm_ν()`](@ref)              | Coefficient ``ν_m`` for the 1-term r.r. boosting ``ℓ → ℓ+1``                   |
| [`CMB.Legendre.Plm_α()`](@ref)              | Coefficient ``α_ℓ^m`` for the 2-term r.r. acting on the ``(ℓ,m)`` term         |
| [`CMB.Legendre.Plm_β()`](@ref)              | Coefficient ``β_ℓ^m`` for the 2-term r.r. acting on the ``(ℓ-1,m)`` term       |

As a concrete example, we'll walk through how ``λ_ℓ^m(x)`` is defined to have the
spherical harmonic normalization baked in.

```math
\begin{align}
    λ_ℓ^m(x) &≡ N_ℓ^m P_ℓ^m(x)
    \\
    N_ℓ^m &= \sqrt{\frac{2ℓ+1}{4π} \frac{(ℓ-m)!}{(ℓ+m)!}}
\end{align}
```

Baking in the normalization happens by changing the coefficients in the recursion
relations given in the [Definitions and Properties](@ref legendre_defn) section.
For our purposes, they take on the form:
```math
\begin{align}
    P_{\ell+1}^m(x) &= \alpha_{\ell+1}^m x P_\ell^m(x)
        - \beta_{\ell+1}^m P_{\ell-1}^m(x)
        \label{eqn:cus_rr_2term}
    \\
    P_{m+1}^{m+1}(x) &= \mu_{m+1} \sqrt{1-x^2} P_m^m(x)
        \label{eqn:cus_rr_1term_lm}
    \\
    P_{m+1}^m(x) &= \nu_m x P_m^m(x)
        \label{eqn:cus_rr_1term_l}
\end{align}
```
The normalization is encoded in the coefficients ``α_ℓ^m``, ``β_ℓ^m``, ``μ_m``, and
``ν_m``.
For the standard (unity) normalization, these take on the values
```math
\begin{align}
    α_ℓ^m &= \frac{2ℓ - 1}{ℓ - m} \\
    β_ℓ^m &= \frac{ℓ + m - 1}{ℓ - m} \\
    μ_m &= 2ℓ - 1 \\
    ν_m &= 2ℓ + 1
\end{align}
```
by simply identifying the coefficients from Eqns.
``\ref{eqn:std_rr_2term}``–``\ref{eqn:std_rr_1term_l}`` on each of the ``P_ℓ^m(x)`` terms
on the right hand side.
For other normalizations, we multiply through by the normalization factor
appropriate for the left-hand side of the equations, rearrange terms to
correctly normalize the terms on the right, and identify the coefficients left
over.
For example, ``α_ℓ^m`` and ``β_ℓ^m`` for ``λ_ℓ^m(x)`` are determined by starting with
Eq. ``\ref{eqn:std_rr_2term}`` and multiply through by ``N_{ℓ+1}^m``.
The left-hand side by definition is ``λ_{ℓ+1}^m``, leaving us with
```math
\begin{align}
    \begin{split}
        λ_{ℓ+1}^m &= \frac{2ℓ + 1}{ℓ - m + 1} x
            \sqrt{\frac{2ℓ+3}{4π} \frac{(ℓ-m+1)!}{(ℓ+m+1)!}} P_ℓ^m(x) -
            \\
            &\quad\quad \frac{ℓ+m}{ℓ-m+1} \sqrt{\frac{2ℓ+3}{4π}
            \frac{(ℓ-m+1)!}{(ℓ+m+1)!}} P_{ℓ-1}^m(x)
    \end{split}
\end{align}
```
Through judicious use of algebra, the terms on the right-hand side can be manipulated
to gather terms of the form ``N_ℓ^m P_ℓ^m(x)`` and ``N_{ℓ-1}^m P_{ℓ-1}^m(x)``, leaving us
with
```math
\begin{align}
    λ_{ℓ+1}^m &= \sqrt{\frac{2ℓ+3}{2ℓ-1} \frac{4ℓ^2 - 1}{(ℓ+1)^2 - m^2}} x
        λ_ℓ^m(x) -
        \sqrt{\frac{2ℓ+3}{2ℓ-1} \frac{ℓ^2 - m^2}{(ℓ+1)^2 - m^2}}
        λ_{ℓ-1}^m(x)
\end{align}
```
We identify each of the two square root terms as ``α_{ℓ+1}^m`` and ``β_{ℓ+1}^m`` since
they are the cofficients appropriate for generating ``λ_{ℓ+1}^m(x)``.
Doing so with the other two recurrence relation equations, we obtain:
```math
\begin{align}
    α_ℓ^m &= \sqrt{\frac{2ℓ+1}{2ℓ-3} \frac{4(ℓ-1)^2 - 1}{ℓ^2 - m^2}} \\
    β_ℓ^m &= \sqrt{\frac{2ℓ+1}{2ℓ-3} \frac{(ℓ-1)^2 - m^2}{ℓ^2 - m^2}} \\
    μ_m &= \sqrt{1 + \frac{1}{2m}} \\
    ν_m &= \sqrt{2m + 3}
\end{align}
```
The final math required is to define the initial condition ``λ_0^0(x)``.
This is straight forward given the definition:
```math
\begin{align}
    λ_0^0(x) &= N_0^0 P_0^0(x) = \sqrt{\frac{1}{4π}} × 1 \\
    λ_0^0(x) &= \sqrt{\frac{1}{4π}}
\end{align}
```

We now have all the information required to define a custom Legendre normalization.
Begin by importing the types and methods which will need to be extended:
```jldoctest λNorm
julia> using CMB.Legendre

julia> import CMB.Legendre: AbstractLegendreNorm, Plm_00, Plm_μ, Plm_ν, Plm_α, Plm_β
```
We'll call our new normalization `λNorm`, which must be a subclass of
`AbstractLegendreNorm`.
```jldoctest λNorm
julia> struct λNorm <: AbstractLegendreNorm end
```
The initial condition is specified by providing a method of `Plm_00` which takes our
normalization trait type as the first argument.
(The second argument can be useful if some extra type information is required to set
up a type-stable algorithm.)
```jldoctest λNorm
julia> Plm_00(::λNorm, T::Type) = sqrt(1 / 4π)
Plm_00 (generic function with 4 methods)
```
Finally, we provide methods which encode the cofficients as well:
```jldoctest λNorm
julia> function Plm_α(::λNorm, T::Type, l::Integer, m::Integer)
           fac1 = (2l + 1) / ((2l - 3) * (l^2 - m^2))
           fac2 = 4*(l-1)^2 - 1
           return sqrt(fac1 * fac2)
       end
Plm_α (generic function with 4 methods)

julia> function Plm_β(::λNorm, T::Type, l::Integer, m::Integer)
           fac1 = (2l + 1) / ((2l - 3) * (l^2 - m^2))
           fac2 = (l-1)^2 - m^2
           return sqrt(fac1 * fac2)
       end
Plm_β (generic function with 4 methods)

julia> Plm_μ(::λNorm, T::Type, m::Integer) = sqrt(1 + 1 / 2m)
Plm_μ (generic function with 4 methods)

julia> Plm_ν(::λNorm, T::Type, m::Integer) = sqrt(3 + 2m)
Plm_ν (generic function with 4 methods)
```

With just those 5 methods provided, the full Legendre framework is available,
including precomputing the coefficients.
```jldoctest λNorm
julia> legendre(λNorm(), 700, 500, 0.4)
0.35366224602811

julia> coeff = LegendreNormCoeff{λNorm,Float64}(700);

julia> legendre(coeff, 700, 500, 0.4)
0.35366224602811
```

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
