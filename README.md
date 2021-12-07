# Spglib

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://singularitti.github.io/Spglib.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://singularitti.github.io/Spglib.jl/dev)
[![Build Status](https://github.com/singularitti/Spglib.jl/workflows/CI/badge.svg)](https://github.com/singularitti/Spglib.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/singularitti/Spglib.jl?svg=true)](https://ci.appveyor.com/project/singularitti/Spglib-jl)
[![Build Status](https://cloud.drone.io/api/badges/singularitti/Spglib.jl/status.svg)](https://cloud.drone.io/singularitti/Spglib.jl)
[![Build Status](https://api.cirrus-ci.com/github/singularitti/Spglib.jl.svg)](https://cirrus-ci.com/github/singularitti/Spglib.jl)
[![pipeline status](https://gitlab.com/singularitti/Spglib.jl/badges/master/pipeline.svg)](https://gitlab.com/singularitti/Spglib.jl/-/pipelines)
[![coverage report](https://gitlab.com/singularitti/Spglib.jl/badges/master/coverage.svg)](https://gitlab.com/singularitti/Spglib.jl/-/jobs)
[![Coverage](https://codecov.io/gh/singularitti/Spglib.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/singularitti/Spglib.jl)
[![Open in Visual Studio Code](https://open.vscode.dev/badges/open-in-vscode.svg)](https://open.vscode.dev/organization/repository)

`Spglib` is a Julia wrapper of the C library [`spglib`](https://github.com/spglib/spglib).
It is used for finding and handling crystal symmetries.
Thanks to Julia's binary artifact mechanism, the installation and usage of it should be
smooth after Julia 1.3.

There was already a package [`LibSymspg.jl`](https://github.com/unkcpz/LibSymspg.jl),
but it is [no longer actively maintained](https://github.com/unkcpz/LibSymspg.jl/issues/4).
And it does not support the latest versions of `spglib`.
It also has some different design decisions with this package, including, but not limited to,
naming conventions, input types, and return types of functions.

The code is [hosted on GitHub](https://github.com/singularitti/Spglib.jl), with some
continuous integration services to test its validity.

## Compatibility

- [Julia version: `v1.3.0` to `v1.7.0`](https://julialang.org/downloads/)
- Dependencies:
  - [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) `v0.8.3` and above
  - [`StructHelpers.jl`](https://github.com/jw3126/StructHelpers.jl) `v0.1.0` and above
  - [`UnPack.jl`](https://github.com/mauro3/UnPack.jl) `v1.0.0` and above
  - [`spglib_jll.jl`](https://github.com/JuliaBinaryWrappers/spglib_jll.jl) `v1.14.1+0` and above
- OS: macOS, Linux, Windows, and FreeBSD
- Architecture: x86, x64, ARM

## Installation

To install `Spglib`, please open Julia's interactive session (known as REPL) and
press `]` key in the REPL to use the [package mode](https://docs.julialang.org/en/v1/stdlib/Pkg/),
then type the following command

For stable release

```julia
(@v1.6) pkg> add Spglib
```

For current master

```julia
(@v1.6) pkg> add Spglib#master
```
