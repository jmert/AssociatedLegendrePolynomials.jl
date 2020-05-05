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

| Interfaces to extend/implement              | Brief description                                                              |
|:------------------------------------------- |:------------------------------------------------------------------------------ |
| [`Legendre.AbstractLegendreNorm`](@ref) | Supertype of normalization trait types                                         |
| [`Legendre.Plm_00()`](@ref)             | Value of ``N_0^0 P_0^0(x)`` for the given normalization                        |
| [`Legendre.Plm_μ()`](@ref)              | Coefficient ``μ_m`` for the 1-term r.r. boosting ``ℓ → ℓ+1`` and ``m → m+1``   |
| [`Legendre.Plm_ν()`](@ref)              | Coefficient ``ν_m`` for the 1-term r.r. boosting ``ℓ → ℓ+1``                   |
| [`Legendre.Plm_α()`](@ref)              | Coefficient ``α_ℓ^m`` for the 2-term r.r. acting on the ``(ℓ,m)`` term         |
| [`Legendre.Plm_β()`](@ref)              | Coefficient ``β_ℓ^m`` for the 2-term r.r. acting on the ``(ℓ-1,m)`` term       |

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
``\ref{eqn:cus_rr_2term}``–``\ref{eqn:cus_rr_1term_l}`` on each of the ``P_ℓ^m(x)`` terms
on the right hand side.
For other normalizations, we multiply through by the normalization factor
appropriate for the left-hand side of the equations, rearrange terms to
correctly normalize the terms on the right, and identify the coefficients left
over.
For example, ``α_ℓ^m`` and ``β_ℓ^m`` for ``λ_ℓ^m(x)`` are determined by starting with
Eq. ``\ref{eqn:cus_rr_2term}`` and multiply through by ``N_{ℓ+1}^m``.
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
julia> using Legendre

julia> import Legendre: AbstractLegendreNorm, Plm_00, Plm_μ, Plm_ν, Plm_α, Plm_β
```
We'll call our new normalization `λNorm`, which must be a subclass of
`AbstractLegendreNorm`.
```jldoctest λNorm
julia> struct λNorm <: AbstractLegendreNorm end
```
The initial condition is specified by providing a method of `Plm_00` which takes our
normalization trait type as the first argument.
(The second argument can be useful if some extra type information is required to set
up a type-stable algorithm, which we'll ignore here for the sake of simplicity.)
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
