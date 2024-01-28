# Spglib

| **Documentation** | [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://singularitti.github.io/Spglib.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://singularitti.github.io/Spglib.jl/dev/)                                                                                                                                                                                                                                                                             |
| :---------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Build Status**  | [![Build Status](https://github.com/singularitti/Spglib.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/singularitti/Spglib.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Build Status](https://ci.appveyor.com/api/projects/status/github/singularitti/Spglib.jl?svg=true)](https://ci.appveyor.com/project/singularitti/Spglib-jl)[![Build Status](https://api.cirrus-ci.com/github/singularitti/Spglib.jl.svg)](https://cirrus-ci.com/github/singularitti/Spglib.jl) |
|   **Coverage**    | [![Coverage](https://github.com/singularitti/Spglib.jl/badges/main/coverage.svg)](https://github.com/singularitti/Spglib.jl/commits/main) [![Coverage](https://codecov.io/gh/singularitti/Spglib.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/singularitti/Spglib.jl)                                                                                                                                                                                                                      |
|    **Others**     | [![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle) [![License](https://img.shields.io/github/license/singularitti/Spglib.jl)](https://github.com/singularitti/Spglib.jl/blob/main/LICENSE)                                                                                                                                                                                                                                       |

The code, which is [hosted on GitHub](https://github.com/singularitti/Spglib.jl), is tested
using various continuous integration services for its validity.

This repository is created and maintained by
[@singularitti](https://github.com/singularitti), and contributions are highly welcome.

## Package features

`Spglib.jl` is a Julia wrapper of the C library [Spglib](https://github.com/spglib/spglib).
See [this page](https://spglib.readthedocs.io/en/latest/interface.html#julia-interface).
It is used for finding and handling crystal symmetries.

There was a package [`LibSymspg.jl`](https://juliahub.com/ui/Packages/LibSymspg/D1i7g)
thanks to [@unkcpz](https://github.com/unkcpz).
However, it is no longer actively maintained.

This package:

- enables the finding and handling of crystal symmetries;
- includes new naming conventions, input types, and return types of functions.

## Installation

The package can be installed with the Julia package manager.
From [the Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/), type `]` to enter
the [Pkg mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode) and run:

```julia-repl
pkg> add Spglib
```

Or, equivalently, via [`Pkg.jl`](https://pkgdocs.julialang.org/v1/):

```julia
julia> import Pkg; Pkg.add("Spglib")
```

## Documentation

- [**STABLE**](https://singularitti.github.io/Spglib.jl/stable/) — **documentation of the most recently tagged version.**
- [**DEV**](https://singularitti.github.io/Spglib.jl/dev/) — _documentation of the in-development version._

## Project status

The package is developed for and tested against Julia `v1.6` and above on Linux, macOS, and
Windows.

## Questions and contributions

You can post usage questions on
[our discussion page](https://github.com/singularitti/Spglib.jl/discussions).

We welcome contributions, feature requests, and suggestions. If you encounter any problems,
please open an [issue](https://github.com/singularitti/Spglib.jl/issues).
The [Contributing](@ref) page has
a few guidelines that should be followed when opening pull requests and contributing code.

## Stargazers over time

[![Stargazers over time](https://starchart.cc/singularitti/Spglib.jl.svg?variant=adaptive)](https://starchart.cc/singularitti/Spglib.jl)
