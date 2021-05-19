```@meta
CurrentModule = Spglib
```

# API

## Types

There are two types, `Dataset` and `SpacegroupType`, correspond to `SpglibDataset`
and `SpglibSpacegroupType`, respectively. They store basic information of a symmetry search.
The struct `Cell` is to contain input data of a symmetry search.

```@docs
Cell
Dataset
SpacegroupType
```

## Methods

Some methods are exported here.
You can find their official documentation [on this page](https://spglib.github.io/spglib/api.html).

```@docs
get_hall_number_from_symmetry
get_dataset
get_spacegroup_type
get_international
get_schoenflies
standardize_cell
find_primitive
refine_cell
niggli_reduce!
delaunay_reduce!
get_multiplicity
get_ir_reciprocal_mesh
get_version
```
