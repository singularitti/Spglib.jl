# Spglib

Documentation for [Spglib](https://github.com/singularitti/Spglib.jl).

See the [Index](@ref main-index) for the complete list of documented functions
and types.

The code, which is [hosted on GitHub](https://github.com/singularitti/Spglib.jl), is tested
using various continuous integration services for its validity.

This repository is created and maintained by
[@singularitti](https://github.com/singularitti), and contributions are highly welcome.

## Package features

- It enables the finding and handling of crystal symmetries.
- The design of the package includes new naming conventions, input types, and return types of functions.

## Installation

The package can be installed with the Julia package manager.
From [the Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/), type `]` to enter
the [Pkg mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode) and run:

```julia-repl
pkg> add Spglib
```

Or, equivalently, via [`Pkg.jl`](https://pkgdocs.julialang.org/v1/):

```@repl
import Pkg; Pkg.add("Spglib")
```

## Documentation

- [**STABLE**](https://singularitti.github.io/Spglib.jl/stable) — **documentation of the most recently tagged version.**
- [**DEV**](https://singularitti.github.io/Spglib.jl/dev) — _documentation of the in-development version._

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

## Manual outline

```@contents
Pages = [
    "man/installation.md",
    "man/definitions.md",
    "man/dataset.md",
    "man/magnetic_dataset.md",
    "man/examples.md",
    "man/troubleshooting.md",
    "developers/contributing.md",
    "developers/style-guide.md",
    "developers/design-principles.md",
]
Depth = 3
```

## Library outline

```@contents
Pages = ["lib/public.md", "lib/internals/error.md", "lib/internals/reciprocal.md"]
```

### [Index](@id main-index)

```@index
Pages = ["lib/public.md"]
```
