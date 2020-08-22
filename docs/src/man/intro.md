# [Introduction](@id intro)

```@contents
Pages = ["intro.md"]
Depth = 2
```

## [Definition and Properties](@id legendre_defn)

The associated Legendre polynomials ``P_\ell^m(x)`` are the solution to the differential
equation
```math
\begin{align}
    (1-x^2) \frac{d^2}{dx^2}P_\ell^m(x) - 2x \frac{d}{dx}P_\ell^m(x) + \left[ \ell(\ell+1) -
        \frac{m^2}{1-x^2} \right] P_\ell^m(x) = 0
\end{align}
```
which arises as the colatitude ``\theta`` part of solving Laplace's equation
``\nabla^2 \psi + \lambda\psi = 0`` in spherical coordinates (where ``x = \cos(\theta)``).

There are several different conventions used to define ``P_\ell^m`` that provide
different properties, but the convention used here is typical of quantum
mechanics and obeys the following properties:

* Solutions only exist for integer ``\ell`` and ``m``, where ``\ell ≤ 0`` and
  ``|m| \le \ell``.

* The associated Legendre functions are normalized such that ``P_0^0`` is unity and have
  orthogonality conditions,
  ```math
  \begin{align}
      \int_{-1}^1 P_\ell^m(x) P_{\ell'}^{m}(x)\,\mathrm{d}x
          = \frac{2}{2\ell+1} \frac{(\ell+m)!}{(\ell-m)!}
          \delta_{\ell\ell'}
  \end{align}
  ```
  for constant ``m`` and
  ```math
  \begin{align}
      \int_{-1}^1 \frac{P_\ell^m(x) P_{\ell}^{m'}(x)}{1-x^2}\,\mathrm{d}x
          = \frac{1}{m} \frac{(\ell+m)!}{(\ell-m)!} \delta_{mm'}
  \end{align}
  ```
  for constant ``\ell``, where ``\delta`` is the Kronecker delta.

* The phase convention for the Legendre functions is chosen such that the negative orders
  are related to positive orders according to,
  ```math
  \begin{align}
      P_\ell^{-m}(x) = (-1)^m \frac{(\ell-m)!}{(\ell+m)!} P_\ell^m(x)
  \end{align}
  ```

* The Legendre functions can be enumerated for non-negative ``m`` using the three
  following recursion relations (given the initial condition ``P_0^0(x)``):
  ```math
  \begin{align}
      P_{\ell+1}^{\ell+1}(x) &= -(2\ell+1)\sqrt{1-x^2} P_\ell^\ell(x)
      \label{eqn:std_rr_1term_lm}
      \\
      P_{\ell+1}^\ell(x) &= x(2\ell+1)P_\ell^\ell(x)
      \label{eqn:std_rr_1term_l}
      \\
      (\ell - m + 1)P_{\ell+1}^m(x) &= (2\ell+1)xP_\ell^m(x) - (\ell+m)P_{\ell-1}^m(x)
      \label{eqn:std_rr_2term}
  \end{align}
  ```


