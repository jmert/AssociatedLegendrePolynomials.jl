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

