export standardize_cell, find_primitive, refine_cell

# See https://github.com/spglib/spglib/blob/444e061/python/spglib/spglib.py#L415-L463 and https://github.com/unkcpz/LibSymspg.jl/blob/f342e72/src/cell-reduce-api.jl#L3-L35
"""
    standardize_cell(cell::AbstractCell, symprec=1e-5; to_primitive=false, no_idealize=false)

Return the standardized cell.

The standardized unit cell (see [Spglib conventions of standardized unit cell](@ref)) is
generated from an input unit cell structure and its symmetry found by
the symmetry search. The choice of the setting for each space group
type is as explained for [`get_dataset`](@ref).
Usually `to_primitive=false` and `no_idealize=false` are recommended to
set and this setting results in the same behavior as
`spg_refine_cell`.

The standardized unit cell is generated from an input unit cell structure and
its symmetry found by the symmetry search. The choice of the setting for each
space group type is as explained for [`get_dataset`](@ref).

# Arguments
- `cell`: the input cell to standardize.
- `symprec`: the tolerance for symmetry search.
- `to_primitive=true` is used to create the standardized primitive cell
  with the transformation matrices shown at
  [Transformation to the primitive cell](@ref), otherwise `to_primitive=false`
  must be specified.
- `no_idealize=false` is used to idealize the lengths and angles of basis
  vectors with adjusting the positions of atoms to nearest exact
  positions according to crystal symmetry. However the crystal can be
  rotated in Cartesian coordinates by the idealization of the basis
  vectors. `no_idealize=true` disables this. The detail of the
  idealization (`no_idealize=false`) is written at
  [Idealization of unit cell structure](@ref idealization). `no_idealize=true` may be useful when we want
  to leave basis vectors and atomic positions in Cartesian coordinates fixed.
"""
function standardize_cell(
    cell::AbstractCell, symprec=1e-5; to_primitive=false, no_idealize=false
)
    lattice, _positions, _atoms = _unwrap_convert(cell)
    num_atom = natoms(cell)
    allocations = 4  # See https://github.com/spglib/spglib/blob/77a8e5d/src/spglib.h#L440
    positions = Matrix{Cdouble}(undef, 3, num_atom * allocations)
    atoms = Vector{Cint}(undef, num_atom * allocations)
    positions[:, 1:num_atom] .= _positions
    atoms[1:num_atom] .= _atoms
    num_atom_std = @ccall libsymspg.spg_standardize_cell(
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        num_atom::Cint,
        to_primitive::Cint,
        no_idealize::Cint,
        symprec::Cdouble,
    )::Cint  # Note: not `num_atom`!
    check_error()
    new_atom_indices = atoms[1:num_atom_std]  # See issue #150
    atoms_record = atomtypes(cell)  # Record the original atoms with labels
    return Cell(
        Lattice(transpose(lattice)),  # We have to `transpose` back because of `_unwrap_convert`!
        collect(eachcol(positions[:, 1:num_atom_std])),
        collect(atoms_record[index] for index in new_atom_indices),
    )
end

"""
    find_primitive(cell::AbstractCell, symprec=1e-5)

Find the primitive cell of an input unit cell.

This function is now a shortcut of `standardize_cell` with `to_primitive=true`
and `no_idealize=false`.
"""
find_primitive(cell::AbstractCell, symprec=1e-5) =
    standardize_cell(cell, symprec; to_primitive=true, no_idealize=false)

"""
    refine_cell(cell::AbstractCell, symprec=1e-5)

Return the refined cell.

This function is now a shortcut of `standardize_cell` with `to_primitive=false`
and `no_idealize=false`.

The standardized crystal structure is obtained from a non-standard crystal
structure which may be slightly distorted within a symmetry recognition
tolerance, or whose primitive vectors are differently chosen, etc.
"""
refine_cell(cell::AbstractCell, symprec=1e-5) =
    standardize_cell(cell, symprec; to_primitive=false, no_idealize=false)
