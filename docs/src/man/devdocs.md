# Developer Documentation

```@meta
DocTestFilters = Regex[
        r"Ptr{0x[0-9a-f]+}",
        r"[0-9\.]+ seconds( \(.*\))?",
        ]
```

```@contents
Pages = ["devdocs.md"]
Depth = 2
```

## [Custom normalizations](@id customnorm)
`Legendre` provides the standard and spherical harmonic normalizations by default, but
arbitrary normalizations are also supported.
The mile-high overview is that the initial condition and recurrence relation (r.r.)
coefficients are all methods which dispatch on a normalization trait type, so a new
normalization is added by simply extending appropriate types and methods.
The following table lists all of the types to extend and method specialization to
implement.

### Normalization Interface

| Interfaces to extend/implement          | Brief description                                                                            |
|:--------------------------------------- |:-------------------------------------------------------------------------------------------- |
| [`Legendre.AbstractLegendreNorm`](@ref) | Supertype of normalization trait types                                                       |
| [`Legendre.initcond()`](@ref)           | Value of ``N_0^0 P_0^0(x)`` for the given normalization                                      |
| [`Legendre.coeff_μ()`](@ref)            | Coefficient ``μ_ℓ`` for the 1-term r.r. boosting ``ℓ-1 → ℓ`` and ``m-1 → m`` where ``m = ℓ`` |
| [`Legendre.coeff_ν()`](@ref)            | Coefficient ``ν_ℓ`` for the 1-term r.r. boosting ``ℓ-1 → ℓ``                                 |
| [`Legendre.coeff_α()`](@ref)            | Coefficient ``α_ℓ^m`` for the 2-term r.r. acting on the ``(ℓ-1,m)`` term                     |
| [`Legendre.coeff_β()`](@ref)            | Coefficient ``β_ℓ^m`` for the 2-term r.r. acting on the ``(ℓ-2,m)`` term                     |

| Optional interfaces                   | Brief description                      |
|:------------------------------------- |:-------------------------------------- |
| [`Legendre.boundscheck_hook()`](@ref) | Hook to participate in bounds checking |


### Example implementation

As a concrete example, we'll walk through how ``λ_ℓ^m(x)`` is defined to have the
spherical harmonic normalization baked in.

```math
\begin{align}
    λ_ℓ^m(x) &≡ N_ℓ^m P_ℓ^m(x)
    \\
    N_ℓ^m &= \sqrt{\frac{2ℓ+1}{4π} \frac{(ℓ-m)!}{(ℓ+m)!}}
\end{align}
```

[^1]:
    Note that here we have shifted the indices by 1 compared to the definitions
    in the introduction such that the left-hand side is always written in terms
    of degree ``ℓ`` rather than ``ℓ+1``.

Baking in the normalization happens by changing the coefficients in the recursion
relations given in the [Definitions and Properties](@ref legendre_defn) section[^1].
For our purposes, they take on the form:
```math
\begin{align}
    P_ℓ^ℓ(x) &= \mu_ℓ \sqrt{1-x^2} P_{ℓ-1}^{ℓ-1}(x)
        \label{eqn:cus_rr_1term_lm}
    \\
    P_ℓ^{ℓ-1}(x) &= \nu_ℓ x P_{ℓ-1}^{ℓ-1}(x)
        \label{eqn:cus_rr_1term_l}
    \\
    P_ℓ^m(x) &= \alpha_ℓ^m x P_{ℓ-1}^m(x)
        - \beta_ℓ^m P_{ℓ-2}^m(x)
        \label{eqn:cus_rr_2term}
\end{align}
```
The normalization is encoded in the coefficients ``μ_ℓ``, ``ν_ℓ``, ``α_ℓ^m``,
``β_ℓ^m``.
For the standard (unity) normalization, these take on the values
```math
\begin{align}
    μ_ℓ &= 2ℓ - 1 \\
    ν_ℓ &= 2ℓ - 1 \\
    α_ℓ^m &= \frac{2ℓ - 1}{ℓ - m} \\
    β_ℓ^m &= \frac{ℓ + m - 1}{ℓ - m}
\end{align}
```
by simply identifying the coefficients from Eqns.
``\ref{eqn:cus_rr_2term}``–``\ref{eqn:cus_rr_1term_l}`` on each of the ``P_ℓ^m(x)`` terms
on the right hand side.

For other normalizations, we multiply through by the normalization factor
appropriate for the left-hand side of the equations, rearrange terms to
correctly normalize the terms on the right, and identify the coefficients left
over.
For example, ``α_ℓ^m`` and ``β_ℓ^m`` for ``λ_ℓ^m(x)`` are determined by starting with
Eq. ``\ref{eqn:cus_rr_2term}`` and multiply through by ``N_ℓ^m``.
The left-hand side by definition is ``λ_ℓ^m``, leaving us with
```math
\begin{align}
    \begin{split}
        λ_ℓ^m &= \frac{2ℓ-1}{ℓ-m} x
            \sqrt{\frac{2ℓ+1}{4π} \frac{(ℓ-m)!}{(ℓ+m)!}} P_{ℓ-1}^m(x) -
            \\
            &\quad\quad \frac{ℓ+m-1}{ℓ-m} \sqrt{\frac{2ℓ+1}{4π}
            \frac{(ℓ-m)!}{(ℓ+m)!}} P_{ℓ-2}^m(x)
    \end{split}
\end{align}
```
Through judicious use of algebra, the terms on the right-hand side can be manipulated
to gather terms of the form ``N_{ℓ-1}^m P_{ℓ-1}^m(x)`` and
``N_{ℓ-2}^m P_{ℓ-2}^m(x)``, leaving us with
```math
\begin{align}
    λ_ℓ^m &= \sqrt{\frac{2ℓ+1}{2ℓ-3} \frac{4(ℓ-1)^2 - 1}{ℓ^2 - m^2}} x
        λ_{ℓ-1}^m(x) -
        \sqrt{\frac{2ℓ+1}{2ℓ-3} \frac{(ℓ-1)^2 - m^2}{ℓ^2 - m^2}}
        λ_{ℓ-2}^m(x)
\end{align}
```
We identify each of the two square root terms as ``α_ℓ^m`` and ``β_ℓ^m`` since
they are the cofficients appropriate for generating ``λ_ℓ^m(x)``.
Doing so with the other two recurrence relation equations, we obtain:
```math
\begin{align}
    μ_ℓ &= \sqrt{1 + \frac{1}{2ℓ}} \\
    ν_ℓ &= \sqrt{2ℓ + 1} \\
    α_ℓ^m &= \sqrt{\frac{2ℓ+1}{2ℓ-3} \frac{4(ℓ-1)^2 - 1}{ℓ^2 - m^2}} \\
    β_ℓ^m &= \sqrt{\frac{2ℓ+1}{2ℓ-3} \frac{(ℓ-1)^2 - m^2}{ℓ^2 - m^2}}
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
julia> using Legendre

julia> import Legendre: AbstractLegendreNorm, initcond, coeff_μ, coeff_ν, coeff_α, coeff_β
```
We'll call our new normalization `λNorm`, which must be a subclass of
`AbstractLegendreNorm`.
```jldoctest λNorm
julia> struct λNorm <: AbstractLegendreNorm end
```
The initial condition is specified by providing a method of `initcond` which takes our
normalization trait type as the first argument.
(The second argument can be useful if some extra type information is required to set
up a type-stable algorithm, which we'll ignore here for the sake of simplicity.)
```jldoctest λNorm
julia> initcond(::λNorm, T::Type) = sqrt(1 / 4π)
initcond (generic function with 4 methods)
```
Finally, we provide methods which encode the cofficients as well:
```jldoctest λNorm
julia> function coeff_α(::λNorm, T::Type, l::Integer, m::Integer)
           fac1 = (2l + 1) / ((2l - 3) * (l^2 - m^2))
           fac2 = 4*(l-1)^2 - 1
           return sqrt(fac1 * fac2)
       end
coeff_α (generic function with 4 methods)

julia> function coeff_β(::λNorm, T::Type, l::Integer, m::Integer)
           fac1 = (2l + 1) / ((2l - 3) * (l^2 - m^2))
           fac2 = (l-1)^2 - m^2
           return sqrt(fac1 * fac2)
       end
coeff_β (generic function with 4 methods)

julia> coeff_μ(::λNorm, T::Type, l::Integer) = sqrt(1 + 1 / 2l)
coeff_μ (generic function with 4 methods)

julia> coeff_ν(::λNorm, T::Type, l::Integer) = sqrt(1 + 2l)
coeff_ν (generic function with 4 methods)
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
