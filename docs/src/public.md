```@meta
CurrentModule = Spglib
```

# API

## Types

There are two types, `Dataset` and `SpacegroupType`, correspond to `SpglibDataset`
and `SpglibSpacegroupType`, respectively. They store basic information of a symmetry search.
The struct `Cell` is to contain input data of a symmetry search.

In `Spglib`, basis vectors are represented by three column vectors, following
the convention of
[`spglib`](https://spglib.github.io/spglib/definition.html#basis-vectors-mathbf-a-mathbf-b-mathbf-c-or-mathbf-a-1-mathbf-a-2-mathbf-a-3).
Coordinates of an atomic point are represented as three fractional values
relative to basis vectors. So when constructing a `Cell`:

- The lattice can be a `3×3` matrix with columns as basis vectors, or it can be a
  vector containing the three basis vectors. To get those basis vectors from a
  `Cell`, use `basis_vectors`.
- The atomic positions can be a `3×N` matrix, where `N` denotes the number of
  atoms in a cell. Or it can be a vector of `N` vectors, where each vector represents an
  atom.
- The `types` variable corresponds to different atomic types.

```@docs
Cell
Dataset
SpacegroupType
```

## Methods

Some methods are exported here.
You can find their official documentation [on this page](https://spglib.github.io/spglib/api.html).

```@docs
basis_vectors
get_symmetry
get_hall_number_from_symmetry
get_dataset
get_dataset_with_hall_number
get_spacegroup_type
get_symmetry_from_database
get_spacegroup_number
get_international
get_schoenflies
standardize_cell
find_primitive
refine_cell
niggli_reduce
delaunay_reduce
get_multiplicity
get_ir_reciprocal_mesh
get_version
```
