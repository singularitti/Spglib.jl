# Definitions and conventions

```@contents
Pages = ["definitions.md"]
Depth = 3
```

!!! warning
Our definitions and conventions are mostly adapted from
[here](https://spglib.readthedocs.io/en/latest/definition.html), with some minor
differences, such as the matrix representation of lattices.

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

in Cartesian coordinates.

Therefore, a lattice is represented as

```math
\mathbf{A} =
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} =
\begin{bmatrix}
    a_x & b_x & c_x \\
    a_y & b_y & c_y \\
    a_z & b_z & c_z
\end{bmatrix}.
```

Depending on the situation,
`\begin{bmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \mathbf{a}_3 \end{bmatrix}`
is used instead of
`\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix}`.

A reciprocal lattice is its inverse, represented as three columns vectors, too:

```math
\mathbf{B} =
\bigl(\mathbf{A}^{-1}\bigr)^\intercal =
\begin{bmatrix} \mathbf{b}_1 & \mathbf{b}_2 & \mathbf{b}_3 \end{bmatrix}
```

so that

```math
\mathbf{A} \mathbf{B}^\intercal = \mathbf{B}^\intercal \mathbf{A} = \mathbf{I},
```

where `\mathbf{I}` is the `3 \times 3` identity matrix.

We choose this convention because it is convenient for converting reduced reciprocal
coordinates `\mathbf{x}^\ast` to Cartesian coordinates using the expression
`\mathbf{B} \mathbf{x}^\ast`.

This is analogous to the convention employed in [Atomic point coordinates](@ref) for the
definition of reduced coordinates.

!!! note
In crystallography, the convention used is
`\mathbf{a}_i \cdot \mathbf{b}_j = \delta_{ij}`, where `\delta_{ij}` is the
[Kronecker delta](https://en.wikipedia.org/wiki/Kronecker_delta).
This is in contrast to the solid-state physics convention, which is
`\mathbf{a}_i \cdot \mathbf{b}_j = 2\pi\delta_{ij}`.

## Atomic point coordinates

Coordinates of an atomic point `\mathbf{x}` are represented
as three fractional values relative to basis vectors as follows,

```math
\mathbf{x} = \begin{bmatrix}
    x_1 \\
    x_2 \\
    x_3
\end{bmatrix},
```

where `0 \le x_i < 1`. A position vector `\mathbf{r}` in
Cartesian coordinates is obtained by

```math
\mathbf{r} = \mathbf{A} \mathbf{x} =
\begin{bmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \mathbf{a}_3 \end{bmatrix}
\begin{bmatrix}
    x_1 \\
    x_2 \\
    x_3
\end{bmatrix}
```

or

```math
\mathbf{r} = \sum_i x_i \mathbf{a}_i.
```

!!! note
In the Python version of Spglib,
lattice parameters `lattice` are given by a `3\times 3` matrix with floating
point values, where `\mathbf{a}`, `\mathbf{b}`, `\mathbf{c}` are
given as rows, which results in the transpose of the definition for
C-API. That is, in Python, the basis vectors are written as
[follows](https://spglib.readthedocs.io/en/latest/variable.html#lattice):

    ```python
    [ [ a_x, b_x, c_x ],
      [ a_y, b_y, c_y ],
      [ a_z, b_z, c_z ] ]
    ```

    Here, we adopt the C-API convention, i.e., writing basis vectors as columns.

## Space group operation and change of basis

### Symmetry operation `(\mathbf{W}, \mathbf{w})`

A symmetry operation consists of a pair of the rotation part
`\mathbf{W}` and translation part `\mathbf{w}`,
and is represented as `(\mathbf{W}, \mathbf{w})`.
The symmetry operation transfers `\mathbf{x}` to
`\tilde{\mathbf{x}}` as follows:

```math
\tilde{\mathbf{x}} = \mathbf{W} \mathbf{x} + \mathbf{w}.
```

### Transformation matrix `\mathbf{P}` and origin shift `\mathbf{p}`

The transformation matrix `\mathbf{P}` changes choice of
basis vectors as follows

```math
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} =
\begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}
\mathbf{P},
```

where `\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix}` and
`\begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}`
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

#### Active/reverse/alibi transformation

A space group operation having no translation part sends an atom to
another point by

```math
\tilde{\mathbf{x}} = \mathbf{W} \mathbf{x},
```

where `\tilde{\mathbf{x}}` and `\mathbf{x}` are
represented with respect to the same basis vectors
`\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix}`.
Equivalently the rotation is achieved by
rotating the basis vectors:

```math
\begin{bmatrix} \tilde{\mathbf{a}} & \tilde{\mathbf{b}} & \tilde{\mathbf{c}} \end{bmatrix} =
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} \mathbf{W}
```

with keeping the representation of the atomic point coordinates
`\mathbf{x}` because

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

## Spglib conventions of standardized unit cell

### Transformation to the primitive cell

In the standardized unit cells, there are five different centring
types available, base centerings of A and C, rhombohedral (R), body-centred (I),
and face-centred (F). The transformation is applied to the
standardized unit cell by

```math
\begin{bmatrix} \mathbf{a}_\text{p} & \mathbf{b}_\text{p} & \mathbf{c}_\text{p} \end{bmatrix} =
\begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}
\mathbf{P}
```

where `\mathbf{a}_\text{p}`, `\mathbf{b}_\text{p}`, and `\mathbf{c}_\text{p}`
are the basis vectors of the primitive cell and `\mathbf{P}` is the
transformation matrix from the standardized unit cell to the primitive
cell. Matrices `\mathbf{P}` for different centring types are given as follows:

```math
\mathbf{P}_\text{A} = \begin{bmatrix}
    1 & 0 & 0 \\
    0 & \dfrac{1}{2} & \dfrac{-1}{2} \\
    0 & \dfrac{1}{2} & \dfrac{1}{2}
\end{bmatrix},
\quad
\mathbf{P}_\text{C} = \begin{bmatrix}
    \dfrac{1}{2} & \dfrac{1}{2} & 0 \\
    \dfrac{-1}{2} & \dfrac{1}{2} & 0 \\
    0 & 0 & 1
\end{bmatrix},
\quad
\mathbf{P}_\text{R} = \begin{bmatrix}
    \dfrac{2}{3} & \dfrac{-1}{3} & \dfrac{-1}{3} \\
    \dfrac{1}{3} & \dfrac{1}{3} & \dfrac{\bar{2}}{3} \\
    \dfrac{1}{3} & \dfrac{1}{3} & \dfrac{1}{3}
\end{bmatrix},
\quad
\mathbf{P}_\text{I} = \begin{bmatrix}
    \dfrac{-1}{2} & \dfrac{1}{2} & \dfrac{1}{2} \\
    \dfrac{1}{2} & \dfrac{-1}{2} & \dfrac{1}{2} \\
    \dfrac{1}{2} & \dfrac{1}{2} & \dfrac{-1}{2}
\end{bmatrix},
\quad
\mathbf{P}_\text{F} = \begin{bmatrix}
    0 & \dfrac{1}{2} & \dfrac{1}{2} \\
    \dfrac{1}{2} & 0 & \dfrac{1}{2} \\
    \dfrac{1}{2} & \dfrac{1}{2} & 0
\end{bmatrix}.
```

The choice of transformation matrix depends on the purpose.

For rhombohedral lattice systems with the H setting (hexagonal lattice),
`\mathbf{P}_\text{R}` is applied to obtain
primitive basis vectors. However, with the R setting (rhombohedral lattice),
no transformation matrix is used because it is already a primitive cell.

### [Idealization of unit cell structure](@id idealization)

Spglib allows tolerance parameters to match a slightly distorted unit
cell structure to a space group type with some higher symmetry. Using
obtained symmetry operations, the distortion is removed to idealize
the unit cell structure. The coordinates of atomic points are
idealized using respective site-symmetries (Grosse-Kunstleve _et al_. (2002)).
The basis vectors are idealized by forcing them into
respective lattice shapes as follows. In this treatment, except for
triclinic crystals, crystals can be rotated in Cartesian coordinates,
which is the different type of transformation from that of the
change-of-basis transformation explained above.

#### Triclinic lattice

- Niggli-reduced cell is used for choosing `\mathbf{a}`, `\mathbf{b}`, and `\mathbf{c}`.
- `\mathbf{a}` is set along `+x` direction of Cartesian coordinates.
- `\mathbf{b}` is set in `x`-`y` plane of Cartesian
  coordinates so that `\mathbf{a}\times\mathbf{b}` is along
  `+z` direction of Cartesian coordinates.

#### Monoclinic lattice

- The `b`-axis is taken as the unique axis.
- `\alpha = 90^\circ` and `\gamma = 90^\circ`, while `90^\circ < \beta < 120^\circ`.
- `\mathbf{a}` is set along `+x` direction of Cartesian coordinates.
- `\mathbf{b}` is set along `+y` direction of Cartesian coordinates.
- `\mathbf{c}` is set in `x`-`z` plane of Cartesian coordinates.

#### Orthorhombic lattice

- `\alpha = \beta = \gamma = 90^\circ`.
- `\mathbf{a}` is set along `+x` direction of Cartesian coordinates.
- `\mathbf{b}` is set along `+y` direction of Cartesian coordinates.
- `\mathbf{c}` is set along `+z` direction of Cartesian coordinates.

#### Tetragonal lattice

- `\alpha = \beta = \gamma = 90^\circ`.
- `a=b`.
- `\mathbf{a}` is set along `+x` direction of Cartesian coordinates.
- `\mathbf{b}` is set along `+y` direction of Cartesian coordinates.
- `\mathbf{c}` is set along `+z` direction of Cartesian coordinates.

#### Rhombohedral lattice

- `\alpha = \beta = \gamma`.
- `a=b=c`.
- Let `\mathbf{a}`, `\mathbf{b}`, and `\mathbf{c}`
  projected on `x`-`y` plane in Cartesian coordinates be
  `\mathbf{a}_{xy}`, `\mathbf{b}_{xy}`, and
  `\mathbf{c}_{xy}`, respectively, and their angles be
  `\alpha_{xy}`, `\beta_{xy}`,
  `\gamma_{xy}`, respectively.
- Let `\mathbf{a}`, `\mathbf{b}`, and `\mathbf{c}`
  projected along `z`-axis in Cartesian coordinates be
  `\mathbf{a}_{z}`, `\mathbf{b}_{z}`, and
  `\mathbf{c}_{z}`, respectively.
- `\mathbf{a}_{xy}` is set along the ray `30^\circ`
  rotated counter-clockwise from the `+x`
  direction of Cartesian coordinates, and `\mathbf{b}_{xy}` and
  `\mathbf{c}_{xy}` are placed by angles `120^\circ` and
  `240^\circ` from `\mathbf{a}_{xy}` counter-clockwise,
  respectively.
- `\alpha_{xy} = \beta_{xy} = \gamma_{xy} = 120^\circ`.
- `a_{xy} = b_{xy} = c_{xy}`.
- `a_{z} = b_{z} = c_{z}`.

#### Hexagonal lattice

- `\alpha = \beta = 90^\circ`, `\gamma = 120^\circ`.
- `a=b`.
- `\mathbf{a}` is set along `+x` direction of Cartesian coordinates.
- `\mathbf{b}` is set in `x`-`y` plane of Cartesian coordinates.
- `\mathbf{c}` is set along `+z` direction of Cartesian coordinates.

#### Cubic lattice

- `\alpha = \beta = \gamma = 90^\circ`.
- `a=b=c`.
- `\mathbf{a}` is set along `+x` direction of Cartesian coordinates.
- `\mathbf{b}` is set along `+y` direction of Cartesian coordinates.
- `\mathbf{c}` is set along `+z` direction of Cartesian coordinates.

### Rotation introduced by idealization

In the idealization step presented above, the input unit cell crystal
structure can be rotated in the Cartesian coordinates. The rotation
matrix `\mathbf{R}` of this rotation is defined by

```math
\begin{bmatrix} \bar{\mathbf{a}}_\text{s} & \bar{\mathbf{b}}_\text{s} & \bar{\mathbf{c}}_\text{s} \end{bmatrix} =
\mathbf{R} \begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}.
```

This rotation matrix rotates the standardized crystal structure before idealization
`\begin{bmatrix} \mathbf{a}_\text{s} & \mathbf{b}_\text{s} & \mathbf{c}_\text{s} \end{bmatrix}`
to that after idealization
`\begin{bmatrix} \bar{\mathbf{a}}_\text{s} & \bar{\mathbf{b}}_\text{s} & \bar{\mathbf{c}}_\text{s} \end{bmatrix}`
in Cartesian coordinates of the given input unit cell.
