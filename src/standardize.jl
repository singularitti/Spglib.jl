# See https://github.com/spglib/spglib/blob/444e061/python/spglib/spglib.py#L415-L463 and https://github.com/unkcpz/LibSymspg.jl/blob/f342e72/src/cell-reduce-api.jl#L3-L35
"""
    standardize_cell(cell::Cell; to_primitive=false, no_idealize=false, symprec=1e-5)

Return standardized cell.

The standardized unit cell is generated from an input unit cell structure and
its symmetry found by the symmetry search. The choice of the setting for each
space group type is as explained for [`get_dataset`](@ref).
"""
function standardize_cell(
    cell::Cell;
    to_primitive = false,
    no_idealize = false,
    symprec = 1e-5,
)
    lattice, positions, types = _expand_cell(cell)
    to_primitive = Base.cconvert(Cint, to_primitive)
    no_idealize = Base.cconvert(Cint, no_idealize)
    num_atom = Base.cconvert(Cint, length(types))
    allocations = 4
    _positions = Matrix{Cdouble}(undef, 3, num_atom * allocations)
    _types = Vector{Cint}(undef, num_atom * allocations)
    _positions[:, 1:num_atom] = positions
    _types[1:num_atom] = types
    num_atom_std = ccall(
        (:spg_standardize_cell, libsymspg),
        Cint,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Cint, Cdouble),
        lattice,
        _positions,
        _types,
        num_atom,
        to_primitive,
        no_idealize,
        symprec,
    )  # Note: not `num_atom`!
    if num_atom_std <= 0
        throw(SpglibError("Cell standardization failed!"))
    end
    # We have to `transpose` back because of `_expand_cell`!
    return Cell(transpose(lattice), collect(eachcol(_positions))[1:num_atom_std], _types[1:num_atom_std])
end

"""
    find_primitive(cell::Cell, symprec=1e-5)

Find the primitive cell of an input unit cell.

This function is now a shortcut of `standardize_cell` with `to_primitive = true`
and `no_idealize = false`.
"""
find_primitive(cell::Cell, symprec = 1e-5) =
    standardize_cell(cell; to_primitive = true, no_idealize = false, symprec = symprec)

"""
    refine_cell(cell::Cell, symprec=1e-5)

Return refined cell.

The standardized crystal structure is obtained from a non-standard crystal
structure which may be slightly distorted within a symmetry recognition
tolerance, or whose primitive vectors are differently chosen, etc.
This function is now a shortcut of `standardize_cell` with `to_primitive = false`
and `no_idealize = false`.
"""
refine_cell(cell::Cell, symprec = 1e-5) =
    standardize_cell(cell; to_primitive = false, no_idealize = false, symprec = symprec)
