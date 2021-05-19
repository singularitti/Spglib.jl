# Installation

```@contents
Pages = ["installation.md"]
Depth = 3
```

To install this package, first, you need to install a `julia` executable from
[its official website](https://julialang.org/downloads/). Version `v1.0.0` and
above is required. This package may not work on `v0.7` and below.

## Installing Julia

### on macOS

If you are using a Mac, and have [Homebrew](https://brew.sh) installed, open
`Terminal.app` and type:

```shell
brew cask install julia
```

### on Linux

On Linux, the best way to install Julia is to use the Generic Linux Binaries. The
[JILL](https://github.com/abelsiqueira/jill) script does this for you. Just run

```shell
bash -ci "$(curl -fsSL https://raw.githubusercontent.com/abelsiqueira/jill/master/jill.sh)"
```

installs Julia into `$HOME/.local/bin`. This script also has a Python version,
[JILL.py](https://github.com/johnnychen94/jill.py). It can also be used on macOS.

## Installing the package

Now I am using [macOS](https://en.wikipedia.org/wiki/MacOS) as a standard
platform to explain the following steps:

1. Open `Terminal.app`, and type `julia` to start a Julia session.

2. Run

   ```julia
   julia> using Pkg; Pkg.update(); Pkg.add("Spglib")
   ```

   and wait for its finish.

3. Run

   ```julia
   julia> using Spglib
   ```

   and have fun!

4. While using, please keep this Julia session alive. Restarting may recompile
   the package and cost some time.

## Reinstall

1. In the same Julia session, run

   ```julia
   julia> Pkg.rm("Spglib"); Pkg.gc()
   ```

2. Press `ctrl+d` to quit the current session. Start a new Julia session and
   repeat the above steps.
