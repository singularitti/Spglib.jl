# Spglib

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://singularitti.github.io/Spglib.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://singularitti.github.io/Spglib.jl/dev)
[![Build Status](https://github.com/singularitti/Spglib.jl/workflows/CI/badge.svg)](https://github.com/singularitti/Spglib.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/singularitti/Spglib.jl?svg=true)](https://ci.appveyor.com/project/singularitti/Spglib-jl)
[![Build Status](https://cloud.drone.io/api/badges/singularitti/Spglib.jl/status.svg)](https://cloud.drone.io/singularitti/Spglib.jl)
[![Build Status](https://api.cirrus-ci.com/github/singularitti/Spglib.jl.svg)](https://cirrus-ci.com/github/singularitti/Spglib.jl)
[![Build Status](https://travis-ci.com/singularitti/Spglib.jl.svg?branch=master)](https://travis-ci.com/singularitti/Spglib.jl)
[![pipeline status](https://gitlab.com/singularitti/Spglib.jl/badges/master/pipeline.svg)](https://gitlab.com/singularitti/Spglib.jl/-/pipelines)
[![Coverage](https://codecov.io/gh/singularitti/Spglib.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/singularitti/Spglib.jl)
[![Coverage Status](https://coveralls.io/repos/github/singularitti/Spglib.jl/badge.svg?branch=master)](https://coveralls.io/github/singularitti/Spglib.jl?branch=master)

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

- [Julia version: `v1.3.0` to `v1.6.1`](https://julialang.org/downloads/)
- Dependencies: see `Project.toml` [`deps` field](Project.toml#L7-L10) and
  [`compat` field](Project.toml#L13-L17)
- OS: macOS, Linux, Windows, and FreeBSD
- Architecture: x86, x64, ARM

## Contributors

This repository is now maintained by [singularitti](https://github.com/singularitti). Thanks to the contribution from [searchengineorientprogramming](https://github.com/searchengineorientprogramming).
