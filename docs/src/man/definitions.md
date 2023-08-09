# Definitions and conventions

```@contents
Pages = ["definitions.md"]
Depth = 3
```

## Basis vectors

In this package, basis vectors are represented by three-column vectors:

```math
\mathbf{a} = \begin{bmatrix}
    a_x \\
    a_y \\
    a_z
\end{bmatrix},
\quad
\mathbf{b} = \begin{bmatrix}
    b_x \\
    b_y \\
    b_z
\end{bmatrix},
\quad
\mathbf{c} = \begin{bmatrix}
    c_x \\
    c_y \\
    c_z
\end{bmatrix},
```

in Cartesian coordinates. Depending on the situation,
``\begin{bmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \mathbf{a}_3 \end{bmatrix}``
is used instead of
``\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix}``.

Therefore, a lattice is represented as

```math
\mathrm{A} =
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} =
\begin{bmatrix}
    a_x & b_x & c_x \\
    a_y & b_y & c_y \\
    a_z & b_z & c_z
\end{bmatrix}.
```

A reciprocal lattice is its inverse, represented as three row vectors:

```math
\mathrm{B} =
\mathrm{A}^{-1} =
\begin{bmatrix}
    \mathbf{b}_1 \\
    \mathbf{b}_2 \\
    \mathbf{b}_3
\end{bmatrix},
```

so that

```math
\mathrm{A} \mathrm{B} = \mathrm{B} \mathrm{A} = \mathrm{I},
```

where ``\mathrm{I}`` is the ``3 \times 3`` identity matrix.

## Crystal coordinates

Coordinates of an atomic point ``\mathbf{x}`` are represented
as three fractional values relative to basis vectors as follows,

```math
\mathbf{x} = \begin{bmatrix}
    x_1 \\
    x_2 \\
    x_3
\end{bmatrix},
```

where ``0 \le x_i < 1``. A position vector ``\mathbf{r}`` in
Cartesian coordinates is obtained by

```math
\mathbf{r} = \mathrm{A} \mathbf{x} = \begin{bmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \mathbf{a}_3 \end{bmatrix} \mathbf{x},
```

or

```math
\mathbf{r} = \sum_i x_i \mathbf{a}_i.
```

## Space group operation and change of basis

### Symmetry operation ``(\mathbf{W}, \mathbf{w})``

A symmetry operation consists of a pair of the rotation part
``\mathbf{W}`` and translation part ``\mathbf{w}``,
and is represented as ``(\mathbf{W}, \mathbf{w})``.
The symmetry operation transfers ``\mathbf{x}`` to
``\tilde{\mathbf{x}}`` as follows:

```math
\tilde{\mathbf{x}} = \mathbf{W} \mathbf{x} + \mathbf{w}.
```

### Transformation matrix ``\mathbf{P}`` and origin shift ``\mathbf{p}``

The transformation matrix ``\mathbf{P}`` changes choice of
basis vectors as follows

```math
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} =
\begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}
\mathbf{P},
```

where ``\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix}`` and
``\begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}``
are the basis vectors of an arbitrary system
and of a standardized system, respectively. In general, the
transformation matrix is not limited for the transformation from the
standardized system, but can be used in between any systems possibly
transformed. It has to be emphasized that the transformation matrix
does not rotate a crystal in Cartesian coordinates, but just
changes the choices of basis vectors.

### Difference between rotation and transformation matrices

A rotation matrix rotates (or mirrors, inverts) the crystal body with
respect to origin. A transformation matrix changes the choice of the
basis vectors, but does not rotate the crystal body.

#### Active/reverse/alibi/ transformation

A space group operation having no translation part sends an atom to
another point by

```math
\tilde{\mathbf{x}} = \mathbf{W} \mathbf{x},
```

where ``\tilde{\mathbf{x}}`` and ``\mathbf{x}`` are
represented with respect to the same basis vectors
``\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix}``.
Equivalently the rotation is achieved by
rotating the basis vectors:

```math
\begin{bmatrix} \tilde{\mathbf{a}} & \tilde{\mathbf{b}} & \tilde{\mathbf{c}} \end{bmatrix} =
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} \mathbf{W} 
```

with keeping the representation of the atomic point coordinates
``\mathbf{x}`` because

```math
\tilde{\mathbf{x}} = \begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix}
\tilde{\mathbf{x}} =
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} \mathbf{W} \mathbf{x} =
\begin{bmatrix} \tilde{\mathbf{a}} & \tilde{\mathbf{b}} & \tilde{\mathbf{c}} \end{bmatrix}
\mathbf{x}.
```

#### Passive/forward/alias transformation

The transformation matrix changes the choice of the basis vectors as:

```math
\begin{bmatrix} \mathbf{a}' & \mathbf{b}' & \mathbf{c}' \end{bmatrix} =
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} \mathbf{P}.
```

The atomic position vector is not altered by this transformation, i.e.,

```math
\begin{bmatrix} \mathbf{a}' & \mathbf{b}' & \mathbf{c}' \end{bmatrix} \mathbf{x}' =
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} \mathbf{x}.
```

However the representation of the atomic point coordinates changes as follows:

```math
\mathbf{P} \mathbf{x}' = \mathbf{x}.
```

because

```math
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} \mathbf{P} \mathbf{x}' =
\begin{bmatrix} \mathbf{a}' & \mathbf{b}' & \mathbf{c}' \end{bmatrix} \mathbf{x}' =
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} \mathbf{x}.
```

## [Transformation to the primitive cell](@id primitive)

In the standardized unit cells, there are five different centring
types available, base centerings of A and C, rhombohedral (R), body-centred (I),
and face-centred (F). The transformation is applied to the
standardized unit cell by

```math
\begin{bmatrix} \mathbf{a}_p & \mathbf{b}_p & \mathbf{c}_p \end{bmatrix} =
\begin{bmatrix} \mathbf{a}_s & \mathbf{b}_s & \mathbf{c}_s \end{bmatrix}
\mathrm{P}
```

where ``\mathbf{a}_p``, ``\mathbf{b}_p``, and ``\mathbf{c}_p``
are the basis vectors of the primitive cell and ``\mathrm{P}`` is the
transformation matrix from the standardized unit cell to the primitive
cell. Matrices ``\mathrm{P}`` for different centring types are given as follows:

```math
\mathrm{P}_\text{A} = \begin{bmatrix}
    1 & 0 & 0 \\
    0 & \dfrac{1}{2} & \dfrac{-1}{2} \\
    0 & \dfrac{1}{2} & \dfrac{1}{2}
\end{bmatrix},
\quad
\mathrm{P}_\text{C} = \begin{bmatrix}
    \dfrac{1}{2} & \dfrac{1}{2} & 0 \\
    \dfrac{-1}{2} & \dfrac{1}{2} & 0 \\
    0 & 0 & 1
\end{bmatrix},
\quad
\mathrm{P}_\text{R} = \begin{bmatrix}
    \dfrac{2}{3} & \dfrac{-1}{3} & \dfrac{-1}{3} \\
    \dfrac{1}{3} & \dfrac{1}{3} & \dfrac{\bar{2}}{3} \\
    \dfrac{1}{3} & \dfrac{1}{3} & \dfrac{1}{3}
\end{bmatrix},
\quad
\mathrm{P}_\text{I} = \begin{bmatrix}
    \dfrac{-1}{2} & \dfrac{1}{2} & \dfrac{1}{2} \\
    \dfrac{1}{2} & \dfrac{-1}{2} & \dfrac{1}{2} \\
    \dfrac{1}{2} & \dfrac{1}{2} & \dfrac{-1}{2}
\end{bmatrix},
\quad
\mathrm{P}_\text{F} = \begin{bmatrix}
    0 & \dfrac{1}{2} & \dfrac{1}{2} \\
    \dfrac{1}{2} & 0 & \dfrac{1}{2} \\
    \dfrac{1}{2} & \dfrac{1}{2} & 0
\end{bmatrix}.
```

The choice of transformation matrix depends on the purpose.

For rhombohedral lattice systems with the H setting (hexagonal lattice),
``\mathrm{P}_\text{R}`` is applied to obtain
primitive basis vectors. However, with the R setting (rhombohedral lattice),
no transformation matrix is used because it is already a primitive cell.
