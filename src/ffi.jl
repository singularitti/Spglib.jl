using UnPack: @unpack

export get_symmetry,
    get_symmetry!,
    get_symmetry_with_collinear_spin,
    get_symmetry_with_collinear_spin!,
    get_hall_number_from_symmetry,
    get_dataset,
    get_spacegroup_type,
    get_international,
    get_schoenflies,
    standardize_cell,
    find_primitive,
    refine_cell,
    niggli_reduce!,
    delaunay_reduce!,
    get_multiplicity,
    get_ir_reciprocal_mesh,
    get_stabilized_reciprocal_mesh,
    get_version

# This is an internal function, do not export!
function get_ccell(cell::Cell{<:AbstractMatrix,<:AbstractMatrix})
    @unpack lattice, positions, types, magmoms = cell
    # Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L32-L35 and https://github.com/spglib/spglib/blob/444e061/python/spglib/spglib.py#L953-L975
    clattice = Base.cconvert(Matrix{Cdouble}, lattice)
    cpositions = Base.cconvert(Matrix{Cdouble}, positions)
    ctypes = Cint[findfirst(isequal(u), unique(types)) for u in types]
    if magmoms !== nothing
        magmoms = Base.cconvert(Vector{Cdouble}, magmoms)
    end
    return Cell(clattice, cpositions, ctypes, magmoms)
end

# Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L70
# This is an internal function, do not export!
function cchars2string(itr)
    vec = collect(Char, Iterators.filter(!iszero, itr))
    return String(vec)
end

# See https://github.com/spglib/spglib/blob/444e061/python/spglib/spglib.py#L115-L165
function get_symmetry(cell::Cell, symprec = 1e-8)
    max_size = length(cell.types) * 48
    rotation = Array{Cint,3}(undef, 3, 3, max_size)
    translation = Array{Cdouble,2}(undef, 3, max_size)
    if cell.magmoms === nothing
        numops = get_symmetry!(rotation, translation, max_size, cell, symprec)
    else
        equivalent_atoms = zeros(length(cell.magmoms))
        primitive_lattice = zeros(Cdouble, 3, 3)
        if ndims(cell.magmoms) == 1
            spin_flips = zeros(length(rotation))
        else
            spin_flips = nothing
        end
        # TODO: unfinished!
    end
    return [
        (rotation = transpose(rotation[:, :, i]), translation = translation[:, i]) for
        i in 1:numops
    ]
end

function get_symmetry!(
    rotation::AbstractArray{T,3},
    translation::AbstractMatrix,
    max_size::Integer,
    cell::Cell,
    symprec = 1e-5,
) where {T}
    @unpack lattice, positions, types = get_ccell(cell)
    rotation = Base.cconvert(Array{Cint,3}, rotation)
    translation = Base.cconvert(Matrix{Cdouble}, translation)
    max_size = Base.cconvert(Cint, max_size)
    number = Base.cconvert(Cint, length(types))
    num_sym = ccall(
        (:spg_get_symmetry, libsymspg),
        Cint,
        (
            Ptr{Cint},
            Ptr{Float64},
            Cint,
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Cint},
            Cint,
            Float64,
        ),
        rotation,
        translation,
        max_size,
        lattice,
        positions,
        types,
        number,
        symprec,
    )
    num_sym == 0 && error("`spg_get_symmetry` failed!")
    return num_sym
end

function get_symmetry_with_collinear_spin!(
    rotation::AbstractArray{T,3},
    translation::AbstractMatrix,
    equivalent_atoms::AbstractVector,
    max_size::Integer,
    cell::Cell,
    symprec = 1e-5,
) where {T}
    @unpack lattice, positions, types, magmoms = get_ccell(cell)
    rotation = Base.cconvert(Array{Cint,3}, rotation)
    translation = Base.cconvert(Matrix{Cdouble}, translation)
    equivalent_atoms = Base.cconvert(Vector{Cint}, equivalent_atoms)
    max_size = Base.cconvert(Cint, max_size)
    number = Base.cconvert(Cint, length(types))
    num_sym = ccall(
        (:spg_get_symmetry_with_collinear_spin, libsymspg),
        Cint,
        (
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cint},
            Cint,
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cint},
            Ptr{Cdouble},
            Cint,
            Cdouble,
        ),
        rotation,
        translation,
        equivalent_atoms,
        max_size,
        lattice,
        positions,
        types,
        magmoms,
        number,
        symprec,
    )
    num_sym == 0 && error("`spg_get_symmetry` failed!")
    return num_sym
end
function get_symmetry_with_collinear_spin(cell::Cell, symprec = 1e-5)
    number = length(cell.types)
    max_size = number * 48
    rotation = Array{Cint,3}(undef, 3, 3, max_size)
    translation = Matrix{Cdouble}(undef, 3, max_size)
    equivalent_atoms = Vector{Cint}(undef, number)
    num_sym = get_symmetry_with_collinear_spin!(
        rotation,
        translation,
        equivalent_atoms,
        max_size,
        cell,
        symprec,
    )
    return rotation[:, :, 1:num_sym], translation[:, 1:num_sym], equivalent_atoms
end

"""
    get_hall_number_from_symmetry(rotation::AbstractArray{T,3}, translation::AbstractMatrix, num_operations::Integer, symprec=1e-5) where {T}

Obtain `hall_number` from the set of symmetry operations.

This is expected to work well for the set of symmetry operations whose
distortion is small. The aim of making this feature is to find space-group-type
for the set of symmetry operations given by the other source than spglib. Note
that the definition of `symprec` is different from usual one, but is given in the
fractional coordinates and so it should be small like `1e-5`.
"""
function get_hall_number_from_symmetry(
    rotation::AbstractArray{T,3},
    translation::AbstractMatrix,
    num_operations::Integer,
    symprec = 1e-5,
) where {T}
    rotation = Base.cconvert(Array{Cint,3}, rotation)
    translation = Base.cconvert(Matrix{Cdouble}, translation)
    num_operations = Base.cconvert(Cint, num_operations)
    return ccall(
        (:spg_get_hall_number_from_symmetry, libsymspg),
        Cint,
        (Ptr{Cint}, Ptr{Float64}, Cint, Float64),
        rotation,
        translation,
        num_operations,
        symprec,
    )
end

"""
    get_multiplicity(cell::Cell, symprec=1e-8)

Return the exact number of symmetry operations. An error is thrown when it fails.
"""
function get_multiplicity(cell::Cell, symprec = 1e-8)
    @unpack lattice, positions, types = get_ccell(cell)
    number = Base.cconvert(Cint, length(types))
    nsymops = ccall(
        (:spg_get_multiplicity, libsymspg),
        Cint,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        lattice,
        positions,
        types,
        number,
        symprec,
    )
    nsymops == 0 && error("Could not determine the multiplicity!")
    return nsymops
end

"""
    get_dataset(cell::Cell, symprec=1e-8)

Search symmetry operations of an input unit cell structure.
"""
function get_dataset(cell::Cell, symprec = 1e-8)
    @unpack lattice, positions, types = get_ccell(cell)
    number = Base.cconvert(Cint, length(types))
    ptr = ccall(
        (:spg_get_dataset, libsymspg),
        Ptr{SpglibDataset},
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        lattice,
        positions,
        types,
        number,
        symprec,
    )
    raw = unsafe_load(ptr)
    return convert(Dataset, raw)
end

"""
    get_spacegroup_type(hall_number)

Translate Hall number to space group type information.
"""
function get_spacegroup_type(hall_number::Integer)
    spgtype = ccall(
        (:spg_get_spacegroup_type, libsymspg),
        SpglibSpacegroupType,
        (Cint,),
        hall_number,
    )
    return convert(SpacegroupType, spgtype)
end

"""
    get_international(cell::Cell, symprec=1e-8)

Return the space group type in Hermann–Mauguin (international) notation.
"""
function get_international(cell::Cell, symprec = 1e-8)
    @unpack lattice, positions, types = get_ccell(cell)
    symbol = Vector{Cchar}(undef, 11)
    exitcode = ccall(
        (:spg_get_international, libsymspg),
        Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        symbol,
        lattice,
        positions,
        types,
        length(types),
        symprec,
    )
    exitcode == 0 && error("Could not determine the international symbol!")
    return cchars2string(symbol)
end

"""
    get_schoenflies(cell::Cell, symprec=1e-8)

Return the space group type in Schoenflies notation.
"""
function get_schoenflies(cell::Cell, symprec = 1e-8)
    @unpack lattice, positions, types = get_ccell(cell)
    symbol = Vector{Cchar}(undef, 7)
    exitcode = ccall(
        (:spg_get_schoenflies, libsymspg),
        Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        symbol,
        lattice,
        positions,
        types,
        length(types),
        symprec,
    )
    exitcode == 0 && error("Could not determine the Schoenflies symbol!")
    return cchars2string(symbol)
end

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
    @unpack lattice, positions, types = get_ccell(cell)
    to_primitive = Base.cconvert(Cint, to_primitive)
    no_idealize = Base.cconvert(Cint, no_idealize)
    number = Base.cconvert(Cint, length(types))
    allocations = 4
    _positions = Matrix{Cdouble}(undef, 3, number * allocations)
    _types = Vector{Cint}(undef, number * allocations)
    _positions[:, 1:number] = positions
    _types[1:number] = types
    num_atom_std = ccall(
        (:spg_standardize_cell, libsymspg),
        Cint,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Cint, Cdouble),
        lattice,
        _positions,
        _types,
        number,
        to_primitive,
        no_idealize,
        symprec,
    )  # Note: not `number`!
    @assert num_atom_std > 0 "standardizing cell failed!"
    return Cell(lattice, _positions[:, 1:num_atom_std], _types[1:num_atom_std])
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

"""
    niggli_reduce!(lattice::AbstractMatrix, symprec=1e-5)

Apply Niggli reduction to input basis vectors `lattice` and the reduced basis vectors are overwritten to `lattice`.
"""
function niggli_reduce!(lattice::AbstractMatrix, symprec = 1e-5)
    clattice = convert(Matrix{Cdouble}, lattice)
    exitcode = ccall(
        (:spg_niggli_reduce, libsymspg),
        Cint,
        (Ptr{Cdouble}, Cdouble),
        clattice,
        symprec,
    )
    iszero(exitcode) && error("Niggli reduce failed!")
    return clattice
end
function niggli_reduce!(cell::Cell, symprec = 1e-5)
    clattice = niggli_reduce!(cell.lattice, symprec)
    cell.lattice[:, :] = clattice
    return cell
end

"""
    delaunay_reduce!(lattice::AbstractMatrix, symprec=1e-5)

Apply Delaunay reduction to input basis vectors `lattice` and the reduced basis vectors are overwritten to `lattice`.
"""
function delaunay_reduce!(lattice::AbstractMatrix, symprec = 1e-5)
    clattice = convert(Matrix{Cdouble}, lattice)
    exitcode = ccall(
        (:spg_delaunay_reduce, libsymspg),
        Cint,
        (Ptr{Cdouble}, Cdouble),
        clattice,
        symprec,
    )
    iszero(exitcode) && error("Delaunay reduce failed!")
    return clattice
end
function delaunay_reduce!(cell::Cell, symprec = 1e-5)
    clattice = delaunay_reduce!(cell.lattice, symprec)
    cell.lattice[:, :] = clattice
    return cell
end

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
returned as the return value.  The time reversal symmetry is
imposed by setting `is_time_reversal`.
"""
function get_ir_reciprocal_mesh(
    cell::Cell,
    mesh,
    is_shift = falses(3);
    is_time_reversal = true,
    symprec = 1e-5,
)
    # Reference: https://github.com/unkcpz/LibSymspg.jl/blob/e912dd3/src/ir-mesh-api.jl#L1-L32
    @assert length(mesh) == length(is_shift) == 3
    @assert all(isone(x) || iszero(x) for x in is_shift)
    # Prepare for input
    @unpack lattice, positions, types = get_ccell(cell)
    mesh = Base.cconvert(Vector{Cint}, mesh)
    is_shift = Base.cconvert(Vector{Cint}, is_shift)
    is_time_reversal = Base.cconvert(Cint, is_time_reversal)
    number = Base.cconvert(Cint, length(types))
    # Prepare for output
    npoints = prod(mesh)
    grid_address = Matrix{Cint}(undef, 3, npoints)  # Julia stores multi-dimensional data in column-major, not row-major (C-style) in memory.
    grid_mapping_table = Vector{Cint}(undef, npoints)
    num_ir = ccall(
        (:spg_get_ir_reciprocal_mesh, libsymspg),
        Cint,
        (
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Cint,
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cint},
            Cint,
            Cdouble,
        ),
        grid_address,
        grid_mapping_table,
        mesh,
        is_shift,
        is_time_reversal,
        lattice,
        positions,
        types,
        number,
        symprec,
    )
    @assert num_ir > 0 "Something wrong happens when finding mesh!"
    return num_ir, grid_mapping_table, grid_address
end

function get_stabilized_reciprocal_mesh(
    rotations::AbstractVector{AbstractMatrix{<:Integer}},
    grid::AbstractVector{<:Integer},
    shift::AbstractVector{<:Integer} = [0, 0, 0];
    qpoints = [[0, 0, 0]],
    is_time_reversal = true,
)
    @assert(length(grid) == length(shift) == 3)
    @assert(all(isone(x) || iszero(x) for x in shift))
    @assert(size(qpoints, 2) == 3)
    npoints = prod(grid)
    grid_address = Matrix{Cint}(undef, npoints, 3)
    mapping = Vector{Cint}(undef, npoints)
    exitcode = ccall(
        (:spg_get_stabilized_reciprocal_mesh, libsymspg),
        Cint,
        (
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Cint,
            Cint,
            Ptr{Cint},
            Cint,
            Ptr{Cint},
        ),
        grid_address,
        mapping,
        grid,
        shift,
        is_time_reversal,
        length(rotations),
        rotations,
        length(qpoints),
        qpoints,
    )
    @assert(exitcode > 0, "Something wrong happens when finding mesh!")
    return mapping, grid_address
end

function Base.convert(::Type{SpacegroupType}, spgtype::SpglibSpacegroupType)
    values = map(fieldnames(SpacegroupType)) do name
        value = getfield(spgtype, name)
        if value isa Cint
            value
        elseif value isa NTuple{N,Cchar} where {N}
            cchars2string(value)
        else  # This should never happen!
            error("unexpected field type $(typeof(value))!")
        end
    end
    return SpacegroupType(values...)
end
