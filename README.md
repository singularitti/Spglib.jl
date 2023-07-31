# Spglib

|                                 **Documentation**                                  |                                                                                                 **Build Status**                                                                                                 |                  **LICENSE**                  |
| :--------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :-------------------------------------------: |
| [![Stable][docs-stable-img]][docs-stable-url] [![Dev][docs-dev-img]][docs-dev-url] | [![Build Status][gha-img]][gha-url] [![Build Status][appveyor-img]][appveyor-url] [![Build Status][cirrus-img]][cirrus-url] [![pipeline status][gitlab-img]][gitlab-url] [![Coverage][codecov-img]][codecov-url] | [![GitHub license][license-img]][license-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://singularitti.github.io/Spglib.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://singularitti.github.io/Spglib.jl/dev
[gha-img]: https://github.com/singularitti/Spglib.jl/workflows/CI/badge.svg
[gha-url]: https://github.com/singularitti/Spglib.jl/actions
[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/singularitti/Spglib.jl?svg=true
[appveyor-url]: https://ci.appveyor.com/project/singularitti/Spglib-jl
[cirrus-img]: https://api.cirrus-ci.com/github/singularitti/Spglib.jl.svg
[cirrus-url]: https://cirrus-ci.com/github/singularitti/Spglib.jl
[gitlab-img]: https://gitlab.com/singularitti/Spglib.jl/badges/master/pipeline.svg
[gitlab-url]: https://gitlab.com/singularitti/Spglib.jl/-/pipelines
[codecov-img]: https://codecov.io/gh/singularitti/Spglib.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/singularitti/Spglib.jl
[license-img]: https://img.shields.io/github/license/singularitti/Spglib.jl
[license-url]: https://github.com/singularitti/Spglib.jl/blob/master/LICENSE

`Spglib` is a Julia wrapper of the C library [`spglib`](https://github.com/spglib/spglib).
It is used for finding and handling crystal symmetries.
Thanks to Julia's binary artifact mechanism, the installation and usage of it should be
smooth after Julia 1.3.

There was already a package [`LibSymspg.jl`](https://juliahub.com/ui/Packages/LibSymspg/D1i7g) by [`@unkcpz`](https://github.com/unkcpz),
but it is [no longer actively maintained](https://github.com/unkcpz/LibSymspg.jl/issues/4).
And it does not support the latest versions of `spglib`.
It also has some different design decisions with this package, including, but not limited to,
naming conventions, input types, and return types of functions.

The code is [hosted on GitHub](https://github.com/singularitti/Spglib.jl), with some
continuous integration services to test its validity.

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add Spglib
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("Spglib")
```

## Compatibility

- [Julia version: `v1.3` to `v1.7`](https://julialang.org/downloads/)
- Dependencies:
  - [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) `v0.8.3` and above
  - [`StructHelpers.jl`](https://github.com/jw3126/StructHelpers.jl) `v0.1.0` and above
  - [`spglib_jll.jl`](https://github.com/JuliaBinaryWrappers/spglib_jll.jl) `v1.15.1+0` and above
- OS: macOS, Linux, Windows, and FreeBSD
- Architecture: x86, x64, ARM

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **documentation of the most recently tagged version.**
- [**DEV**][docs-dev-url] &mdash; _documentation of the in-development version._

## Project Status

The package is tested against, and being developed for, Julia `1.6` and above on Linux,
macOS, and Windows.

## Questions and Contributions

Usage questions can be posted on [our discussion page][discussions-url].

Contributions are very welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems. The [contributing](@ref) page has
a few guidelines that should be followed when opening pull requests and contributing code.

[discussions-url]: https://github.com/singularitti/Spglib.jl/discussions
[issues-url]: https://github.com/singularitti/Spglib.jl/issues
[contrib-url]: https://github.com/singularitti/Spglib.jl/discussions
