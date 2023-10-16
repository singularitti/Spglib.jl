# Examples

```@contents
Pages = ["examples.md"]
Depth = 2
```

## Creating a `Cell`

To create a `Cell`, we first need to create a `Lattice`.
Then we can add atoms and their positions (in crystal coordinates):

```@repl cell
using Spglib
lattice = [
    -3.0179389205999998 -3.0179389205999998 0.0000000000000000
    -5.2272235447000002 5.2272235447000002 0.0000000000000000
    0.0000000000000000 0.0000000000000000 -9.7736219469000005
]
positions = [[2 / 3, 1 / 3, 1 / 4], [1 / 3, 2 / 3, 3 / 4]]
atoms = [1, 1]
cell = Cell(lattice, positions, atoms)
```

## Computing rigid rotation introduced by idealization

This example is from
[here](https://spglib.readthedocs.io/en/latest/definition.html#computing-rigid-rotation-introduced-by-idealization).

In this package, rigid rotation is purposely introduced in the idealization step though this
is unlikely as a crystallographic operation.

```@repl std
using StaticArrays, Spglib
lattice = Lattice([
    [5.0759761474456697, 5.0759761474456697, 0],
    [-2.8280307701821314, 2.8280307701821314, 0],
    [0, 0, 8.57154746],
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
atoms = fill(35, length(positions));
cell = Cell(lattice, positions, atoms)
dataset = get_dataset(cell, 1e-5)
dataset.international_symbol
dataset.spacegroup_number
dataset.transformation_matrix
```

We can see the transformation matrix from the given lattice to the standardized lattice
is the identity matrix, i.e., the given lattice is already a standardized lattice.

```@repl std
std_lattice_before_idealization = convert(Matrix{Float64}, lattice) * inv(dataset.transformation_matrix)
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

Here, ``\mathbf{P}`` is the `dataset.transformation_matrix`.

Note that in contrast to the Python code:

```python
std_lattice_before_idealization = np.dot(
    np.transpose(lattice),
    np.linalg.inv(dataset['transformation_matrix'])).T
```

where there are multiple transpose operations, we do not have to do that in our Julia
code since we choose a column-major order of stacking lattice vectors as described in
[Basis vectors](@ref), and we return transformation matrix in column-major order, too.

Now, we obtain the standardized basis vectors after idealization
``\begin{bmatrix} \bar{\mathbf{a}}_\text{s} & \bar{\mathbf{b}}_\text{s} & \bar{\mathbf{c}}_\text{s} \end{bmatrix}``:

```@repl std
std_lattice_after_idealization = dataset.std_lattice
```

This is different from the standardized basis vectors before idealization
``\begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}``.
Unless this crystal structure is distorted from the crystal structure that has the ideal
symmetry, this means that the crystal was rotated rigidly in the idealization step by

```math
\begin{bmatrix} \bar{\mathbf{a}}_\text{s} & \bar{\mathbf{b}}_\text{s} & \bar{\mathbf{c}}_\text{s} \end{bmatrix} =
\mathbf{R} \begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix},
```

as stated in [Rotation introduced by idealization](@ref).
where ``\mathbf{R}`` is the rotation matrix. This is computed by

```math
\mathbf{R} =
\begin{bmatrix} \bar{\mathbf{a}}_\text{s} & \bar{\mathbf{b}}_\text{s} & \bar{\mathbf{c}}_\text{s} \end{bmatrix}
\bigl(\begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}\bigr)^{-1}
```

In Julia code, this is

```@repl std
𝐀 = convert(Matrix{Float64}, std_lattice_after_idealization)
𝐁 = convert(Matrix{Float64}, std_lattice_before_idealization)
𝐑 = 𝐀 * inv(𝐁)
```

Note also the transpose is not applied here in contrast to the Python code:

```python
R = np.dot(dataset['std_lattice'].T, np.linalg.inv(std_lattice_before_idealization.T))
```

This equals to

```math
\begin{bmatrix}
    \cos \theta & -\sin \theta & 0\\
    \sin \theta & \cos \theta & 0\\
    0 & 0 & 1
\end{bmatrix}
```

where ``\theta = -\pi/4`` and ``\det(\mathbf{R}) = 1`` when no distortion:

```@repl std
θ = -π/4
[
    cos(θ) -sin(θ) 0
    sin(θ) cos(θ) 0
    0 0 1
]
```

Compared to `dataset.std_rotation_matrix`:

```@repl std
dataset.std_rotation_matrix
```

we have approximately the same result.

In summary of the two steps,

```math
\begin{align}
    \begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix} =
    \begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} \mathbf{P}^{-1} &=
    \mathbf{R}^{-1}
    \begin{bmatrix} \bar{\mathbf{a}}_\text{s} & \bar{\mathbf{b}}_\text{s} & \bar{\mathbf{c}}_\text{s} \end{bmatrix},\\
    \begin{bmatrix} \bar{\mathbf{a}}_\text{s} & \bar{\mathbf{b}}_\text{s} & \bar{\mathbf{c}}_\text{s} \end{bmatrix}
    \mathbf{P} &=
    \mathbf{R} \begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}.
\end{align}
```
