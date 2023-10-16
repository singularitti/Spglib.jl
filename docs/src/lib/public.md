# Public API

## Contents

```@contents
Pages = ["public.md"]
Depth = 2
```

## Index

```@index
Pages = ["public.md"]
```

## Public interface

To utilize the public interface, first import the required packages:

```@repl
using CrystallographyCore, Spglib
```

For documentation on `CrystallographyCore.jl`, refer to
[this link](https://mineralscloud.github.io/CrystallographyCore.jl/stable/).

For extended functionalities, consider exploring my other packages:
[`CrystallographyBase.jl`](https://github.com/MineralsCloud/CrystallographyBase.jl),
[`Crystallography.jl`](https://github.com/MineralsCloud/Crystallography.jl),
and [`MillerIndices.jl`](https://github.com/MineralsCloud/MillerIndices.jl).

### Types

```@autodocs
Modules = [Spglib]
Private = false
Order = [:type]
```

### Functions

You can find their official documentation [on this page](https://spglib.github.io/spglib/api.html).

```@autodocs
Modules = [Spglib]
Private = false
Order = [:function]
```
