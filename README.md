# AssociatedLegendrePolynomials.jl — Calculating Associated Legendre Polynomials

| **Documentation**                                                         | **Build Status**                                             |
|:-------------------------------------------------------------------------:|:------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][ci-img]][ci-url][![][codecov-img]][codecov-url] |

AssociatedLegendrePolynomials.jl is a library for computing the Associated Legendre Polynomials.

Design goals of this package include:

  * Native Julia implementation of core routines.

  * Numerical stability and efficiency.

  * Parallelism and efficient memory sharing.

### Installation and usage

Installation and loading is as easy as:

```
pkg> add AssociatedLegendrePolynomials

julia> using AssociatedLegendrePolynomials

# or on julia >= v1.6, importing to a shorter name is possible:

julia> import AssociatedLegendrePolynomials as Legendre

julia> using .Legendre
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://jmert.github.io/AssociatedLegendrePolynomials.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://jmert.github.io/AssociatedLegendrePolynomials.jl/dev

[ci-img]: https://github.com/jmert/AssociatedLegendrePolynomials.jl/actions
[ci-url]: https://github.com/jmert/AssociatedLegendrePolynomials.jl/workflows/CI/badge.svg

[codecov-img]: https://codecov.io/gh/jmert/AssociatedLegendrePolynomials.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/jmert/AssociatedLegendrePolynomials.jl
