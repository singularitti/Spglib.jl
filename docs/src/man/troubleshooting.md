# Troubleshooting

```@contents
Pages = ["troubleshooting.md"]
Depth = 2
```

This page collects some possible errors you may encounter along with tips on how to fix them.
If you have some questions about how to use this code, you are welcome to
[discuss with us](https://github.com/singularitti/Spglib.jl/discussions).

If you have additional tips, please either
[report an issue](https://github.com/singularitti/Spglib.jl/issues/new) or
[submit a pull request](https://github.com/singularitti/Spglib.jl/compare) with suggestions.

## Cannot find the Julia executable

Make sure you have Julia installed in your environment. Please download the latest
[stable version](https://julialang.org/downloads/#current_stable_release) for your platform.
If you are using a *nix system, the recommended way is to use
[Juliaup](https://github.com/JuliaLang/juliaup). If you do not want to install Juliaup
or you are using other platforms that Julia supports, download the corresponding binaries.
Then, create a symbolic link to the Julia executable.
If the path is not in your [`$PATH` environment variable](https://en.wikipedia.org/wiki/PATH_(variable)),
export it to your `$PATH`.

Some clusters, like
[Comet](https://www.sdsc.edu/support/user_guides/comet.html),
or [Expanse](https://www.sdsc.edu/services/hpc/expanse/index.html),
already have Julia installed as a module, you may
just `module load julia` to use it. If not, either install by yourself or contact your
administrator.

See [Installation Guide](@ref) for more information.

## Julia starts slow

First, we recommend you download the latest version of Julia. Usually, the newest version
has the best performance.

If you need to use Julia for a simple, one-time task, you can start the Julia REPL with

```bash
julia --compile=min
```

to minimize compilation or

```bash
julia --optimize=0
```

to minimize optimizations, or just use both. Or you could make a system image
and run with

```bash
julia --sysimage custom-image.so
```

See [Fredrik Ekre's talk](https://youtu.be/IuwxE3m0_QQ?t=313) for details.

## Returned cell symmetry is wrong

Check whether you set the lattice correctly. This is the part where errors can easily occur,
as we adopt a different convention from the Python and C versions.

For example, [the example](https://github.com/spglib/spglib/blob/v2.1.0-rc2/README.md)
shown in Spglib's official documentation is written as follows:

```c
#include <assert.h>
#include "spglib.h"

int main(void) {
    SpglibDataset *dataset;
    // Wurtzite structure (P6_3mc)
    double lattice[3][3] = {
        {3.111, -1.5555, 0}, {0, 2.6942050311733885, 0}, {0, 0, 4.988}};
    double position[4][3] = {
        {1.0 / 3, 2.0 / 3, 0.0},
        {2.0 / 3, 1.0 / 3, 0.5},
        {1.0 / 3, 2.0 / 3, 0.6181},
        {2.0 / 3, 1.0 / 3, 0.1181},
    };
    int types[4] = {1, 1, 2, 2};
    int num_atom = 4;
    double symprec = 1e-5;
    dataset = spg_get_dataset(lattice, position, types, num_atom, symprec);
    assert(dataset->spacegroup_number == 186);
    spg_free_dataset(dataset);
}
```

Thus, the Python correspondence of the code should be:

```python
import numpy as np
import spglib

lattice = np.array([
    [3.111, 0, 0],
    [-1.5555, 2.6942050311733885, 0],
    [0, 0, 4.988]
])
positions = [
    [1.0 / 3, 2.0 / 3, 0.0],
    [2.0 / 3, 1.0 / 3, 0.5],
    [1.0 / 3, 2.0 / 3, 0.6181],
    [2.0 / 3, 1.0 / 3, 0.1181]
]
types = [1, 1, 2, 2]
num_atom = 4
symprec = 1e-5
cell = (lattice, positions, types)
dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)
assert dataset['number'] == 186
```

Note that in Python, the `lattice` is transposed, as explained in its
[official documentation](https://spglib.readthedocs.io/en/latest/variable.html#lattice).

However, the corresponding code in Julia should be written as follows:

```@repl example
using Spglib

lattice = [[3.111, 0, 0], [-1.5555, 2.6942050311733885, 0], [0, 0, 4.988]];
positions = [
    [1.0 / 3, 2.0 / 3, 0.0],
    [2.0 / 3, 1.0 / 3, 0.5],
    [1.0 / 3, 2.0 / 3, 0.6181],
    [2.0 / 3, 1.0 / 3, 0.1181],
];
atoms = [1, 1, 2, 2];
cell = Cell(lattice, positions, atoms)
dataset = get_dataset(cell, 1e-5)
dataset.spacegroup_number
```

Although the Julia definition of our [Basis vectors](@ref) is not transposed
(like the C-API), when written one by one, it still resembles a transposed
version of the lattice in C. However, they both represent the following matrix:

```math
\begin{bmatrix}
    3.111 & -1.5555 & 0 \\
    0 & 2.6942050311733885 & 0 \\
    0 & 0 & 4.988
\end{bmatrix}
```

Of course, you can construct the lattice directly using its matrix form:

```@repl example
lattice = [
    3.111  -1.5555  0.0
    0.0  2.6942050311733885  0.0
    0.0  0.0  4.988
];
```

If you are not careful when writing these matrices, you may encounter unexpected results:

```python
import spglib

lattice = [
    [3.111, -1.5555, 0],
    [0, 2.6942050311733885, 0],
    [0, 0, 4.988]
]
positions = [
    [1.0 / 3, 2.0 / 3, 0.0],
    [2.0 / 3, 1.0 / 3, 0.5],
    [1.0 / 3, 2.0 / 3, 0.6181],
    [2.0 / 3, 1.0 / 3, 0.1181]
]
types = [1, 1, 2, 2]
num_atom = 4
symprec = 1e-5
cell = (lattice, positions, types)
dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)

>>> dataset['number']
4
>>> dataset['international']
'P2_1'
```

Or this:

```@repl example
lattice = [[3.111, -1.5555, 0], [0, 2.6942050311733885, 0], [0, 0, 4.988]];
cell = Cell(lattice, positions, atoms)
dataset = get_dataset(cell, 1e-5)
dataset.spacegroup_number
dataset.international_symbol
```

## Julia's results are different from Python's

For the same reason, the returned results from Python, especially lattices, are
transposed versions of those in Julia.

```python
>>> dataset['primitive_lattice']
array([[ 3.111     ,  0.        ,  0.        ],
       [-1.5555    ,  2.69420503,  0.        ],
       [ 0.        ,  0.        ,  4.988     ]])
>>> dataset['std_lattice']
array([[ 3.111     ,  0.        ,  0.        ],
       [-1.5555    ,  2.69420503,  0.        ],
       [ 0.        ,  0.        ,  4.988     ]])
>>> std_lattice_before_idealization = np.dot(
        np.transpose(lattice),
        np.linalg.inv(dataset['transformation_matrix'])
    ).T
array([[ 3.111     ,  0.        ,  0.        ],
       [-1.5555    ,  2.69420503,  0.        ],
       [ 0.        ,  0.        ,  4.988     ]])
```

These are equivalent to the following Julia matrices:

```@repl
[
    3.111   0.0      0.0
    -1.5555  2.69421  0.0
    0.0     0.0      4.988
];
[
    3.111   0.0      0.0
    -1.5555  2.69421  0.0
    0.0     0.0      4.988
];
[
    3.111   0.0      0.0
    -1.5555  2.69421  0.0
    0.0     0.0      4.988
];
```

However, the actual Julia results should be:

```@repl example
using Spglib: Lattice
dataset.primitive_lattice
dataset.std_lattice
std_lattice_before_idealization =
    convert(Matrix{Float64}, Lattice(cell)) * inv(dataset.transformation_matrix)
```
