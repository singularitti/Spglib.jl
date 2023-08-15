export get_ir_reciprocal_mesh, get_stabilized_reciprocal_mesh

# Doc from https://github.com/spglib/spglib/blob/d1cb3bd/src/spglib.h#L424-L439
"""
    get_ir_reciprocal_mesh(cell::Cell, mesh, is_shift=falses(3); is_time_reversal=true, symprec=1e-5)

Return k-points mesh and k-point map to the irreducible k-points.

Irreducible reciprocal grid points are searched from uniform
mesh grid points specified by `mesh` and `is_shift`.
`mesh` stores three integers. Reciprocal primitive vectors
are divided by the number stored in `mesh` with (0,0,0) point
centering. The centering can be shifted only half of one mesh
by setting `1` or `true` for each `is_shift` element. If `0` or `false` is set for
`is_shift`, it means there is no shift. This limitation of
shifting enables the irreducible k-point search significantly
faster when the mesh is very dense.

The reducible uniform grid points are returned in reduced
coordinates as `grid_address`. A map between reducible and
irreducible points are returned as `grid_mapping_table` as in the indices of
`grid_address`. The number of the irreducible k-points are
returned as the return value. The time reversal symmetry is
imposed by setting `is_time_reversal`.

!!! compat "Version 0.2"
    The returned mapping table is indexed starting at `1`, not `0` as in Python or C.
"""
function get_ir_reciprocal_mesh(
    cell::AbstractCell, mesh, symprec=1e-5; is_shift=falses(3), is_time_reversal=true
)
    # Reference: https://github.com/unkcpz/LibSymspg.jl/blob/e912dd3/src/ir-mesh-api.jl#L1-L32
    if !(length(mesh) == length(is_shift) == 3)
        throw(DimensionMismatch("`mesh` & `is_shift` must be both length-three vectors!"))
    end
    @assert all(isone(x) || iszero(x) for x in is_shift)
    # Prepare for input
    lattice, positions, atoms = _expand_cell(cell)
    mesh = Base.cconvert(Vector{Cint}, mesh)
    is_shift = Base.cconvert(Vector{Cint}, is_shift)
    # Prepare for output
    num_k = prod(mesh)
    grid_address = Matrix{Cint}(undef, 3, num_k)  # Julia stores multi-dimensional data in column-major, not row-major (C-style) in memory.
    ir_mapping_table = Vector{Cint}(undef, num_k)
    num_ir = @ccall libsymspg.spg_get_ir_reciprocal_mesh(
        grid_address::Ptr{Cint},
        ir_mapping_table::Ptr{Cint},
        mesh::Ptr{Cint},
        is_shift::Ptr{Cint},
        is_time_reversal::Cint,
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        natoms(cell)::Cint,
        symprec::Cdouble,
    )::Cint
    check_error()
    ir_mapping_table .+= 1  # See https://github.com/singularitti/Spglib.jl/issues/56
    return num_ir, ir_mapping_table, grid_address
end

function get_stabilized_reciprocal_mesh(
    rotations, mesh, qpoints=[[0, 0, 0]], is_shift=falses(3); is_time_reversal=true
)
    if !(length(mesh) == length(is_shift) == 3)
        throw(DimensionMismatch("`grid` & `is_shift` must be both length-three vectors!"))
    end
    @assert all(isone(x) || iszero(x) for x in is_shift)
    @assert size(qpoints, 2) == 3
    nk = prod(mesh)
    grid_address = Matrix{Cint}(undef, 3, nk)
    ir_mapping_table = Vector{Cint}(undef, nk)
    qpoints = map(Base.Fix1(convert, Vector{Cint}), qpoints)
    nir = @ccall libsymspg.spg_get_stabilized_reciprocal_mesh(
        grid_address::Ptr{Cint},
        ir_mapping_table::Ptr{Cint},
        mesh::Ptr{Cint},
        is_shift::Ptr{Cint},
        is_time_reversal::Cint,
        length(rotations)::Cint,
        rotations::Ptr{Cint},
        length(qpoints)::Cint,
        qpoints::Ptr{Cdouble},
    )::Cint
    check_error()
    ir_mapping_table .+= 1  # See https://github.com/singularitti/Spglib.jl/issues/56
    return nir, ir_mapping_table, grid_address
end
