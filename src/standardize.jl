# See https://github.com/spglib/spglib/blob/444e061/python/spglib/spglib.py#L415-L463 and https://github.com/unkcpz/LibSymspg.jl/blob/f342e72/src/cell-reduce-api.jl#L3-L35
"""
    standardize_cell(cell::Cell; to_primitive=false, no_idealize=false, symprec=1e-5)

Return standardized cell.

The standardized unit cell is generated from an input unit cell structure and
its symmetry found by the symmetry search. The choice of the setting for each
space group type is as explained for [`get_dataset`](@ref).
"""
function standardize_cell(
    cell::AbstractCell; to_primitive=false, no_idealize=false, symprec=1e-5
)
    lattice, _positions, _atoms = _expand_cell(cell)
    n = natoms(cell)
    allocations = 4
    positions = Matrix{Cdouble}(undef, 3, n * allocations)
    atoms = Vector{Cint}(undef, n * allocations)
    positions[:, 1:n] .= _positions
    atoms[1:n] .= _atoms
    n_std = @ccall libsymspg.spg_standardize_cell(
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        n::Cint,
        to_primitive::Cint,
        no_idealize::Cint,
        symprec::Cdouble,
    )::Cint  # Note: not `n`!
    check_error()
    # We have to `transpose` back because of `_expand_cell`!
    return Cell(
        Lattice(transpose(lattice)), collect(eachcol(positions[:, 1:n_std])), atoms[1:n_std]
    )
end

"""
    find_primitive(cell::Cell, symprec=1e-5)

Find the primitive cell of an input unit cell.

This function is now a shortcut of `standardize_cell` with `to_primitive = true`
and `no_idealize = false`.
"""
find_primitive(cell::AbstractCell, symprec=1e-5) =
    standardize_cell(cell; to_primitive=true, no_idealize=false, symprec=symprec)

"""
    refine_cell(cell::Cell, symprec=1e-5)

Return refined cell.

The standardized crystal structure is obtained from a non-standard crystal
structure which may be slightly distorted within a symmetry recognition
tolerance, or whose primitive vectors are differently chosen, etc.
This function is now a shortcut of `standardize_cell` with `to_primitive = false`
and `no_idealize = false`.
"""
refine_cell(cell::AbstractCell, symprec=1e-5) =
    standardize_cell(cell; to_primitive=false, no_idealize=false, symprec=symprec)
