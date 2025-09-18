# Style Guide

```@contents
Pages = ["design-principles.md"]
Depth = 2:3
```

This section describes the coding style rules that apply to our code and that
we recommend you to use it also.

In some cases, our style guide diverges from Julia's official
[Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/) (Please read it!).
All such cases will be explicitly noted and justified.

Our style guide adopts many recommendations from the
[BlueStyle](https://github.com/invenia/BlueStyle).
Please read the [BlueStyle](https://github.com/invenia/BlueStyle)
before contributing to this package.
If these guidelines are not followed, your pull requests may not be accepted.

!!! info
    The style guide is always a work in progress, and not all Spglib code
    follows the rules. When modifying Spglib, please fix the style violations
    of the surrounding code (i.e., leave the code tidier than when you
    started). If large changes are needed, consider separating them into
    another pull request.

## Formatting

### Run JuliaFormatter

Spglib uses [JuliaFormatter](https://github.com/domluna/JuliaFormatter.jl) as
an auto-formatting tool.

We use the options contained in [`.JuliaFormatter.toml`](https://github.com/singularitti/Spglib.jl/blob/main/.JuliaFormatter.toml).

To format your code, `cd` to the Spglib directory, then run:

```julia-repl
julia> using Pkg

julia> Pkg.add("JuliaFormatter")

julia> using JuliaFormatter: format

julia> format("docs"); format("src"); format("test")
```

!!! info
    A continuous integration check verifies that all PRs made to Spglib have
    passed the formatter.

The following sections outline extra style guide points that are not fixed
automatically by JuliaFormatter.

### Use the Julia extension for Visual Studio Code

Please use [Visual Studio Code](https://code.visualstudio.com/) with the
[Julia extension](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia)
to edit, format, and test your code.
For the time being, we do not recommend using editors other than Visual Studio Code to edit your code.

This extension already has [JuliaFormatter](https://github.com/domluna/JuliaFormatter.jl)
integrated. So to format your code, follow the steps listed
[here](https://www.julia-vscode.org/docs/stable/userguide/formatter/).
