export get_ir_reciprocal_mesh, get_stabilized_reciprocal_mesh, eachpoint

"""
    BrillouinZoneMesh(mesh, is_shift, ir_mapping_table, grid_address)

Represent the Brillouin zone mesh for searching (irreducible) reciprocal grid points from
uniform mesh grid points specified by `mesh` and `is_shift`.

Reciprocal primitive vectors are divided by the number stored in `mesh` with ``(0, 0,
0)``-centering. The center of grid mesh is shifted half of a grid spacing along the
corresponding reciprocal axis by setting `1` or `true` to each `is_shift` element. If `0` or
`false` is set to each `is_shift` element, there is no shift. This limitation of shifting
enables the irreducible k-point search to be significantly faster when the mesh is very
dense.

The reducible uniform grid points are stored in fractional coordinates as `grid_address`. A
map between reducible and irreducible points is stored as `ir_mapping_table` in the indices
of `grid_address`.

!!! warning
    The mapping table's indices start from `1`, not `0` as in Python or C.
"""
struct BrillouinZoneMesh
    mesh::SVector{3,Int64}
    is_shift::SVector{3,Bool}
    ir_mapping_table::Vector{Int64}
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

# Arguments
- `cell`: the input cell.
- `mesh`: the mesh numbers along each reciprocal axis. It is given by three integers.
- `symprec`: the tolerance for symmetry search.
- `is_shift`: a 3-Boolean vector. When `is_shift` is set for each reciprocal primitive axis,
  the mesh is shifted along the axis in half of adjacent mesh points irrespective of the
  mesh numbers.
- `is_time_reversal`: whether to impose the time reversal symmetry.

See also [`BrillouinZoneMesh`](@ref).
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
    lattice, positions, atoms = _unwrap_convert(cell)
    mesh = Base.cconvert(Vector{Cint}, mesh)  # Very important to convert!
    is_shift = Base.cconvert(Vector{Cint}, is_shift)  # Very important to convert!
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
    return BrillouinZoneMesh(
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
This function can be used to obtain all mesh grid points by setting
`rotations = [[1 0 0; 0 1 0; 0 0 1]]`, and `qpoints = [[0, 0, 0]]`.

See also [`BrillouinZoneMesh`](@ref).
"""
function get_stabilized_reciprocal_mesh(
    rotations, mesh, qpoints=[[0, 0, 0]]; is_shift=falses(3), is_time_reversal=true
)
    if !(length(mesh) == length(is_shift) == 3)
        throw(DimensionMismatch("`grid` & `is_shift` must be both length-three vectors!"))
    end
    @assert all(isone(x) || iszero(x) for x in is_shift)
    rotations = Base.cconvert(Array{Cint,3}, cat(transpose.(rotations)...; dims=3))
    mesh = Base.cconvert(Vector{Cint}, mesh)  # Very important to convert!
    is_shift = Base.cconvert(Vector{Cint}, is_shift)  # Very important to convert!
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
    return BrillouinZoneMesh(
        SVector{3}(mesh),
        SVector{3}(is_shift),
        ir_mapping_table,
        map(SVector{3,Int64}, eachcol(grid_address)),
    )
end

"""
    eachpoint(result::BrillouinZoneMesh, ir_only=true)

Iterate over the points in the Brillouin zone mesh, with the option to include only
irreducible k-points or all k-points.

See also [`BrillouinZoneMesh`](@ref).

!!! note
    This function only returns an iterator, not a vector of points. To get a vector, use
    `collect(eachpoint(result, ir_only))`.
"""
function eachpoint(result::BrillouinZoneMesh, ir_only=true)
    mesh, shift, grid_address = result.mesh, result.is_shift .// 2, result.grid_address  # true / 2 = 0.5, false / 2 = 0
    if ir_only  # Return only irreducible k-points
        return Iterators.map(unique(result.ir_mapping_table)) do i
            (grid_address[i] .+ shift) .// mesh
        end
    else  # Return all k-points
        return Iterators.map(grid_address) do point
            (point .+ shift) .// mesh
        end
    end
end
