# Magnetic dataset (experimental)

```@contents
Pages = ["definitions.md"]
Depth = 2
```

```@setup dataset
using Spglib
```

The dataset is accessible through the C-structure given by

```@example dataset
dump(MagneticDataset)
```

## Magnetic space-group type

### `uni_number`

The serial number from ``1`` to ``1651`` of UNI or BNS symbols.

### `msg_type`

Magnetic space groups (MSG) is classified by its family space group (FSG) and maximal space
subgroup (XSG). FSG is a non-magnetic space group obtained by ignoring time-reversal term in
MSG. XSG is a space group obtained by picking out non time-reversal operations in MSG.

- `msg_type=1` (type-I): MSG, XSG, FSG are all isomorphic.
- `msg_type=2` (type-II): XSG and FSG are isomorphic, and MSG is generated from XSG and
  pure time reversal operations
- `msg_type=3` (type-III): XSG is a proper subgroup of MSG with isomorphic translational
  subgroups.
- `msg_type=4` (type-IV): XSG is a proper subgroup of MSG with isomorphic point group.

### `hall_number`



### `tensor_rank`

0 for collinear spins, 1 for non-collinear spins

## Magnetic symmetry operations

### `n_operations`

Number of magnetic symmetry operations.

### `rotations`

Rotation (matrix) parts of symmetry operations

### `translations`


Translation (vector) parts of symmetry operations

### `time_reversals`

Time reversal part of magnetic symmetry operations. 1 indicates time reversal operation, and
0 indicates an ordinary operation.

## Symmetrically equivalent atoms

### `n_atoms` and `equivalent_atoms`



## Transformation to standardized setting

### `transformation_matrix` and `origin_shift`

See [Transformation matrix and origin shift](@ref).

## Standardized magnetic crystal structure after idealization

### `n_std_atoms`, `std_lattice`, `std_types`, `std_positions`, and `std_rotation_matrix`

See [Standardized crystal structure after idealization](@ref).

### `std_tensors`



## Intermediate data in symmetry search

### `primitive_lattice`


