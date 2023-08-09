# Examples

```@contents
Pages = ["examples.md"]
Depth = 2
```

## Creating a `Cell`

To create a `Cell`, we first need to create a `Lattice`.
There are multiple ways of doing it. For example, if we know the six lattice constants,
we can do

```@repl 1
using Spglib, Unitful, UnitfulAtomic
lattice₁ = Lattice(4u"nm", 180u"bohr", 3u"angstrom", 90, 90, 90)
```

Or, equivalently,

```@repl 1
lattice₁ == Lattice([
    4u"nm" 0u"m" 0.0u"cm"
    0u"cm" 180.0u"bohr" 0u"m"
    0u"bohr" 0u"nm" (3//1)*u"angstrom"
])
```

Then we can add atoms and their positions (in crystal coordinates):

```@repl 1
lattice₂ = [
    -3.0179389205999998 -3.0179389205999998 0.0000000000000000
    -5.2272235447000002 5.2272235447000002 0.0000000000000000
    0.0000000000000000 0.0000000000000000 -9.7736219469000005
]
positions = [[2 / 3, 1 / 3, 1 / 4], [1 / 3, 2 / 3, 3 / 4]]
atoms = [1, 1]
cell = Cell(lattice₂, positions, atoms)
```

## Computing rigid rotation introduced by idealization

This example is from
[here](https://spglib.readthedocs.io/en/latest/definition.html#computing-rigid-rotation-introduced-by-idealization).

In this package, rigid rotation is purposely introduced in the idealization step though this
is unlikely as a crystallographic operation.

```@repl 1
lattice = Lattice([
    [5.0759761474456697, 5.0759761474456697, 0],  # a
    [-2.8280307701821314, 2.8280307701821314, 0],  # b
    [0, 0, 8.57154746],  # c
]);
positions = [
    [0.0, 0.84688439, 0.1203133],
    [0.0, 0.65311561, 0.6203133],
    [0.0, 0.34688439, 0.3796867],
    [0.0, 0.15311561, 0.8796867],
    [0.5, 0.34688439, 0.1203133],
    [0.5, 0.15311561, 0.6203133],
    [0.5, 0.84688439, 0.3796867],
    [0.5, 0.65311561, 0.8796867],
];
numbers = fill(35, length(positions))
cell = Cell(lattice, positions, numbers)
dataset = get_dataset(cell, 1e-5)
dataset.international_symbol
dataset.spacegroup_number
dataset.transformation_matrix
std_lattice_after_idealization = dataset.std_lattice
std_lattice_before_idealization = Matrix(lattice) * inv(dataset.transformation_matrix)
```

This is based on formula in
[Transformation matrix ``\mathbf{P}`` and origin shift ``\mathbf{p}``](@ref)
and [Passive/forward/alias transformation](@ref):

```math
\begin{align}
    \begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} &=
    \begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}
    \mathbf{P},\\
    \therefore \begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix} &=
    \begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} \mathbf{P}^{-1}.
\end{align}
```

Note that in contrast to the Python code:

```python
std_lattice_before_idealization = np.dot(
    np.transpose(lattice),
    np.linalg.inv(dataset['transformation_matrix'])).T
```

where there are multiple transpose operations, we do not have to do that in our Julia
code since we choose a column-major order of stacking lattice vectors as described in
[Basis vectors](@ref), and we return transformation matrix in column-major order, too.
