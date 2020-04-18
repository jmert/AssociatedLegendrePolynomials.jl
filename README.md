# Legendre.jl — Calculating Associated Legendre Polynomials

| **Documentation**                                                         | **Build Status**                                             |
|:-------------------------------------------------------------------------:|:------------------------------------------------------------:|
| <!--[![][docs-stable-img]][docs-stable-url]--> [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url][![][codecov-img]][codecov-url] |

Legendre.jl is a library for computing the Associated Legendre Polynomials.

Design goals of this package include:

  * Native Julia implementation of core routines.

  * Numerical stability and efficiency.

  * Parallelism and efficient memory sharing.

### Installation and usage

This library is **not** registered in Julia's [General registry][General.jl],
so the package must be installed either by cloning it directly:

```
(@v1.4) pkg> add https://github.com/jmert/Legendre.jl
```

or by making use of my [personal registry][Registry.jl]:

```
(@v1.4) pkg> registry add https://github.com/jmert/Registry.jl
(@v1.4) pkg> add Legendre
```

After installing, just load like any other Julia package:

```
julia> using Legendre
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://jmert.github.io/Legendre.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://jmert.github.io/Legendre.jl/dev

[travis-img]: https://travis-ci.com/jmert/Legendre.jl.svg?branch=master
[travis-url]: https://travis-ci.com/jmert/Legendre.jl

[codecov-img]: https://codecov.io/gh/jmert/Legendre.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/jmert/Legendre.jl

[General.jl]: https://github.com/JuliaRegistries/General
[Registry.jl]: https://github.com/jmert/Registry.jl
