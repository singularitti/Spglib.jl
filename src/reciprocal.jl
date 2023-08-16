export get_ir_reciprocal_mesh, get_stabilized_reciprocal_mesh, eachpoint

struct MeshResult
    mesh::SVector{3,UInt64}
    is_shift::SVector{3,Bool}
    ir_mapping_table::Vector{UInt64}
    grid_address::Vector{SVector{3,Int64}}
end

# Doc from https://github.com/spglib/spglib/blob/d1cb3bd/src/spglib.h#L424-L439
"""
    get_ir_reciprocal_mesh(cell::AbstractCell, mesh, symprec=1e-5; is_shift=falses(3), is_time_reversal=true)

Search irreducible reciprocal grid points from uniform mesh grid points specified by `mesh` and `is_shift`.

Reciprocal primitive vectors are divided by the number stored in `mesh` with ``(0, 0, 0)``-centering.
The center of grid mesh is shifted half of a grid spacing along corresponding reciprocal axis
by setting `1` or `true` to each `is_shift` element. If `0` or `false` is set to each
`is_shift` element, there is no shift. This limitation of
shifting enables the irreducible k-point search significantly
faster when the mesh is very dense.

The reducible uniform grid points are returned in fractional
coordinates as `grid_address`. A map between reducible and
irreducible points are returned as `ir_mapping_table` as in the indices of
`grid_address`. The number of the irreducible k-points are also
returned. The time reversal symmetry is imposed by setting `is_time_reversal=true`.

# Arguments
- `cell`: the input cell.
- `mesh`: the mesh numbers along each reciprocal axis. It is given by three integers.
- `symprec`: the tolerance for symmetry search.
- `is_shift`: a 3-Boolean vector. When `is_shift` is set for each reciprocal primitive axis,
  the mesh is shifted along the axis in half of adjacent mesh points irrespective of the
  mesh numbers.
- `is_time_reversal`: whether to impose the time reversal symmetry.

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
    @ccall libsymspg.spg_get_ir_reciprocal_mesh(
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
    return MeshResult(
        SVector{3}(mesh),
        SVector{3}(is_shift),
        ir_mapping_table,
        map(SVector{3,Int64}, eachcol(grid_address)),
    )
end

"""
    get_stabilized_reciprocal_mesh(rotations, mesh, qpoints=[[0, 0, 0]]; is_shift=falses(3), is_time_reversal=true)

Search irreducible k-points from unique k-point mesh grids
from direct (real space) basis vectors and a set of rotation parts of
symmetry operations in direct space with one or multiple
stabilizers.

The stabilizers are written in fractional coordinates.
Symmetrically equivalent k-points
(stars) in fractional coordinates are stored in `map` as indices of
`grid_address`. The number of reduced k-points with the stabilizers
are returned as the return value.

This function can be used to obtain all mesh grid points by setting
`rotations = [[1 0 0; 0 1 0; 0 0 1]]`, and `qpoints = [[0, 0, 0]]`.
"""
function get_stabilized_reciprocal_mesh(
    rotations, mesh, qpoints=[[0, 0, 0]]; is_shift=falses(3), is_time_reversal=true
)
    if !(length(mesh) == length(is_shift) == 3)
        throw(DimensionMismatch("`grid` & `is_shift` must be both length-three vectors!"))
    end
    @assert all(isone(x) || iszero(x) for x in is_shift)
    rotations = Base.cconvert(Array{Cint,3}, cat(transpose.(rotations)...; dims=3))
    num_k = prod(mesh)
    grid_address = Matrix{Cint}(undef, 3, num_k)
    ir_mapping_table = Vector{Cint}(undef, num_k)
    qpoints = Base.cconvert(Matrix{Cdouble}, reduce(hcat, qpoints))
    @ccall libsymspg.spg_get_stabilized_reciprocal_mesh(
        grid_address::Ptr{Cint},
        ir_mapping_table::Ptr{Cint},
        mesh::Ptr{Cint},
        is_shift::Ptr{Cint},
        is_time_reversal::Cint,
        size(rotations, 3)::Cint,  # `num_rot` in the C-API
        rotations::Ptr{Cint},
        size(qpoints, 2)::Cint,  # `num_q` in the C-API
        qpoints::Ptr{Cdouble},
    )::Cint
    check_error()
    ir_mapping_table .+= 1  # See https://github.com/singularitti/Spglib.jl/issues/56
    return MeshResult(
        SVector{3}(mesh),
        SVector{3}(is_shift),
        ir_mapping_table,
        map(SVector{3,Int64}, eachcol(grid_address)),
    )
end

function eachpoint(result::MeshResult, ir_only=true)
    mesh, shift, grid_address = result.mesh, result.is_shift ./ 2, result.grid_address  # true / 2 = 0.5, false / 2 = 0
    if ir_only
        return Iterators.map(unique(result.ir_mapping_table)) do i
            (grid_address[i] .+ shift) ./ mesh
        end
    else
        return Iterators.map(grid_address) do point
            (point .+ shift) ./ mesh  # Add 1 because `mapping` index starts from 0
        end
    end
end
