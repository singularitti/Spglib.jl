```@meta
CurrentModule = Spglib
```

# Spglib

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

## Manual Outline

```@contents
Pages = [
    "installation.md",
    "portability.md",
    "api.md",
]
Depth = 3
```

## Index

```@index
```
