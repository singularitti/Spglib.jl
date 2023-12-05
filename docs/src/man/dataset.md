# Spglib dataset

```@contents
Pages = ["definitions.md"]
Depth = 2
```

```@setup dataset
using Spglib
```

The dataset is accessible through the `struct` given by

```@example dataset
dump(Dataset)
```

## Space group type

### `spacegroup_number`

The space group type number defined in International Tables for
Crystallography (ITA).

### `hall_number`

The serial number from ``1`` to ``530`` which are found at
[list of space groups (Seto's web site)](https://yseto.net/?page_id=29%3E). Be
sure that this is not a standard crystallographic definition as far as
the author of Spglib knows.

### `international_symbol`

The (full) Hermann–Mauguin notation of space group type is given by .

### `hall_symbol`

The Hall symbol is stored here.

### `choice`

The information on unique axis, setting or cell choices.

## Symmetry operations

### `rotations`, `translations`, and `n_operations`

The symmetry operations of the input unit cell are stored in
`rotations` and `translations`. A crystallographic symmetry
operation ``(\mathbf{W}, \mathbf{w})`` is made from a pair
of rotation ``\mathbf{W}`` and translation
``\mathbf{w}`` parts with the same index. Number of symmetry
operations is given as `n_operations`. The detailed explanation of
the values is found at [`get_symmetry`](@ref).

## Wyckoff positions and symmetrically equivalent atoms

### `n_atoms`

Number of atoms in the input unit cell. This gives the numbers of
elements in `wyckoffs` and `equivalent_atoms`.

### `wyckoffs`

This gives the information of Wyckoff letters by integer
numbers, where ``0``, ``1``, ``2``, ``\ldots``, represent the Wyckoff letters
of ``a``, ``b``, ``c``, ``\ldots`` These are assigned to all atomic positions
of the input unit cell in this order. Therefore the number of elements in
`wyckoffs` is same as the number of atoms in the input unit cell,
which is given by `n_atoms`.

This is determined from the symmetry of the primitive cell.

### `site_symmetry_symbols`

This gives site-symmetry symbols. These are valid for the standard
settings. For different settings and choices belonging to the same
space group type, the same set of the symbols is returned.

This is determined from the symmetry of the primitive cell.

### `equivalent_atoms`

This gives the mapping table from the atomic indices of the input unit
cell to the atomic indices of symmetrically independent atom, such as
`[1, 1, 1, 1, 5, 5, 5, 5]`, where the symmetrically independent
atomic indices are ``1`` and
``5``. We can see that the atoms from ``1`` to ``4`` are mapped to ``1`` and those
from ``5`` to ``8`` are mapped to ``5``.  The number of elements in
`equivalent_atoms` is same as the number of atoms in the input unit
cell, which is given by `n_atoms`.

!!! warning
    You may notice that the indices here differ from those in
    [Spglib's official documentation](https://spglib.readthedocs.io/en/latest/dataset.html#equivalent-atoms),
    where the indices start from ``0``. This discrepancy arises because indices in Julia
    start from ``1`` by default. Consequently, all indices here are incremented by ``1``.

Symmetry operations found for the input cell are used to determine the
equivalent atoms. `equivalent_atoms` and `crystallographic_orbits`
are almost equivalent, but they can be different in a special
case as written in [`get_symmetry`](@ref).

### `crystallographic_orbits`

This is almost equivalent to `equivalent_atoms`. But symmetry of the
primitive cell is used to determine the symmetrically equivalent atoms.

## Transformation matrix and origin shift

### `transformation_matrix` and `origin_shift`

`transformation_matrix` (``\mathbf{P}``) and
`origin_shift` (``\mathbf{p}``) are obtained as a result of
space-group-type matching under a set of unique axis, setting and cell
choices. These are operated to the basis vectors and atomic point
coordinates of the input unit cell as

```math
\begin{align}
    \begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix} &=
    \begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} \mathbf{P}^{-1},\\
    \mathbf{x}_\text{s} &= \mathbf{P} \mathbf{x} + \mathbf{p} \ (\mathrm{mod}\ \mathbf{1}),
\end{align}
```

by which the basis vectors are transformed to those of a
standardized unit cell. Atomic point coordinates are shifted so that
symmetry operations have one of possible standard origins. The
detailed definition is presented at
[Definitions and conventions](@ref).

## Standardized crystal structure after idealization

### `n_std_atoms`, `std_lattice`, `std_types`, and `std_positions`

The standardized crystal structure after [idealization](@ref idealization)
corresponding to a Hall symbol is stored in `n_std_atoms`, `std_lattice`, `std_types`, and
`std_positions`. These output usually contains the rotation in Cartesian coordinates and
rearrangement of the order atoms with respect to the input unit cell.

### `std_rotation_matrix`

Rotation matrix that rotates the standardized crystal structure
before idealization
``\begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}``
to that after idealization
``\begin{bmatrix} \bar{\mathbf{a}}_\text{s} & \bar{\mathbf{b}}_\text{s} & \bar{\mathbf{c}}_\text{s} \end{bmatrix}``
in Cartesian coordinates of the given input unit cell. The rotation
matrix ``\mathbf{R}`` is defined by

```math
\begin{bmatrix} \bar{\mathbf{a}}_\text{s} & \bar{\mathbf{b}}_\text{s} & \bar{\mathbf{c}}_\text{s} \end{bmatrix} =
\begin{bmatrix} \mathbf{R} \mathbf{a}_\text{s} & \mathbf{R} \mathbf{b}_\text{s} & \mathbf{R} \mathbf{c}_\text{s} \end{bmatrix}
```

More precisely, this rotation matrix is an orthonormal matrix. Since
``\begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}``
can be distored, to make ``\mathbf{R}`` orthonormal, it is calculated as
follows. Make cubes of
``\begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}``
and
``\begin{bmatrix} \bar{\mathbf{a}}_\text{s} & \bar{\mathbf{b}}_\text{s} & \bar{\mathbf{c}}_\text{s} \end{bmatrix}``
by

```math
\mathbf{L} = \begin{bmatrix}
\dfrac{\mathbf{a}}{\lvert\mathbf{a}\rvert} &
\dfrac{(\mathbf{a} \times \mathbf{b}) \times \mathbf{a}}{\lvert(\mathbf{a} \times \mathbf{b}) \times \mathbf{a}\rvert} &
\dfrac{\mathbf{a} \times \mathbf{b}}{\lvert\mathbf{a} \times \mathbf{b}\rvert}
\end{bmatrix}.
```

Watching ``\mathbf{L}_\text{s}`` as ``3\times 3`` matrices, ``\mathbf{R}`` is
obtained by solving

```math
\bar{\mathbf{L}}_\text{s} = \mathbf{R} \mathbf{L}_\text{s}.
```

### `std_mapping_to_primitive`

This gives a list of
atomic indices in the primitive cell of the standardized crystal
structure, where the same number presents the same atom in the
primitive cell. By collective the atoms having the same number, a set
of relative lattice points in the standardized crystal structure
is obtained.

## Crystallographic point group

### `pointgroup_symbol`

`pointgroup_symbol` is the symbol of the crystallographic point
group in the Hermann–Mauguin notation. There are 32 crystallographic
point groups

```math
\{1,\ \bar{1},\ 2,\ m,\ 2/m,\ 222,\ mm2,\ mmm,\ 4,\ \bar{4},\ 4/m,\ 422,\ 4mm,\ \bar{4}2m,\ 4/mmm,\ 3,\ \bar{3},\ 32,\ 3m,\ \bar{3}m,\ 6,\ \bar{6},\ 6/m,\ 622,\ 6mm,\ \bar{6}m2,\ 6/mmm,\ 23,\ m\bar{3},\ 432,\ \bar{4}3m,\ m\bar{3}m\}
```

## Intermediate data in symmetry search

A primitive cell is searched from the translational symmetry. This
primitive cell is given by `primitive_lattice` and
`mapping_to_primitive` below.

### `primitive_lattice`

Non-standardized basis vectors of a primitive cell in the input
cell.

### `mapping_to_primitive`

This
gives a list of atomic indices in the primitive cell of the input
crystal structure, where the same number presents the same atom in the
primitive cell. By collective the atoms having the same number, a set
of relative lattice points in the input crystal structure is obtained.
