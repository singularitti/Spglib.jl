# Public API

## Contents

```@contents
Pages = ["public.md"]
Depth = 2
```

## Index

```@index
Pages = ["public.md"]
```

## Public interface

### Types

```@docs
Lattice
SpglibCell
Dataset
SpacegroupType
```

### Functions

You can find their official documentation [on this page](https://spglib.github.io/spglib/api.html).

```@docs
basis_vectors
natoms
atomtypes
get_symmetry
get_symmetry_from_database
get_dataset
get_dataset_with_hall_number
get_multiplicity
get_international
get_schoenflies
get_spacegroup_type
get_spacegroup_type_from_symmetry
get_hall_number_from_symmetry
standardize_cell
find_primitive
refine_cell
niggli_reduce
delaunay_reduce
get_ir_reciprocal_mesh
get_stabilized_reciprocal_mesh
get_version
```
