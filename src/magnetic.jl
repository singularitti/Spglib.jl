export MagneticDataset,
    MagneticSpacegroupType,
    get_symmetry_with_collinear_spin,
    get_symmetry_with_site_tensors,
    get_magnetic_dataset,
    get_magnetic_symmetry_from_database,
    get_magnetic_spacegroup_type,
    get_magnetic_spacegroup_type_from_symmetry,
    is_spin_collinear

"""
    is_spin_collinear(cell::SpglibCell)

Check if the spins (magmoms) in the given `cell` are collinear.
"""
function is_spin_collinear(cell::SpglibCell)
    if isempty(cell.magmoms)
        throw(ArgumentError("`cell.magmoms` is empty!"))
    elseif eltype(cell.magmoms) <: AbstractVector &&
        all(length(site) == 3 for site in cell.magmoms)
        return false
    else
        return true
    end
end

# Python version: https://github.com/spglib/spglib/blob/42527b0/python/spglib/spglib.py#L182-L319
"""
    get_symmetry_with_collinear_spin(cell::SpglibCell, symprec=1e-5)

Find symmetry operations with collinear polarizations (spins) on atoms.

Except for the magmoms in the `cell`, the usage is basically the same as [`get_symmetry`](@ref).
But as an output, `equivalent_atoms` are obtained as the last returned value. The size of this
array is the same of number of atoms in the cell.
"""
function get_symmetry_with_collinear_spin(cell::SpglibCell, symprec=1e-5)
    lattice, positions, atoms, magmoms = _unwrap_convert(cell)
    num_atom = length(cell.magmoms)
    # See https://github.com/spglib/spglib/blob/42527b0/python/spglib/spglib.py#L270
    max_size = 96num_atom  # 96 = 48 × 2 since we have spins
    rotations = Array{Cint,3}(undef, 3, 3, max_size)
    translations = Matrix{Cdouble}(undef, 3, max_size)
    equivalent_atoms = Vector{Cint}(undef, num_atom)
    num_sym = @ccall libsymspg.spg_get_symmetry_with_collinear_spin(
        rotations::Ptr{Cint},
        translations::Ptr{Cdouble},
        equivalent_atoms::Ptr{Cint},
        max_size::Cint,
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        magmoms::Ptr{Cdouble},
        num_atom::Cint,
        symprec::Cdouble,
    )::Cint
    check_error()
    rotations = map(
        SMatrix{3,3,Int32,9} ∘ transpose, eachslice(rotations[:, :, begin:num_sym]; dims=3)
    )  # Remember to transpose, see https://github.com/singularitti/Spglib.jl/blob/8aed6e0/src/core.jl#L195-L198
    translations = map(SVector{3,Float64}, eachcol(translations[:, begin:num_sym]))
    return rotations, translations, equivalent_atoms .+ 1
end

"""
    get_symmetry_with_site_tensors(cell::SpglibCell, symprec=1e-5; with_time_reversal=true)

Return magnetic symmetry operations represented by `rotation`, `translation`, and `spin_flips`.

Returned `spin_flips` represents sign of site tensors after applying time-reversal operations:
``1`` for non time reversal, and ``-1`` for time reversal.
"""
function get_symmetry_with_site_tensors(
    cell::SpglibCell, symprec=1e-5; with_time_reversal=true
)
    lattice, positions, atoms, magmoms = _unwrap_convert(cell)
    num_atom = natoms(cell)
    # See https://github.com/spglib/spglib/blob/42527b0/python/spglib/spglib.py#L270
    max_size = 96num_atom  # 96 = 48 × 2 since we have spins
    rotations = Array{Cint,3}(undef, 3, 3, max_size)
    translations = Matrix{Cdouble}(undef, 3, max_size)
    equivalent_atoms = Vector{Cint}(undef, num_atom)
    primitive_lattice = zeros(Cdouble, 3, 3)
    spin_flips = zeros(Cint, max_size)
    tensor_rank = ndims(magmoms) - 1  # See https://github.com/spglib/spglib/blob/v2.1.0/python/spglib/spglib.py#L275-L276 & https://github.com/spglib/spglib/blob/v2.1.0/python/spglib/spglib.py#L615
    is_axial = iszero(tensor_rank) ? false : true  # Collinear spin & non-collinear spin
    num_sym = @ccall libsymspg.spg_get_symmetry_with_site_tensors(
        rotations::Ptr{Cint},
        translations::Ptr{Cdouble},
        equivalent_atoms::Ptr{Cint},
        primitive_lattice::Ptr{Cdouble},
        spin_flips::Ptr{Cint},
        max_size::Cint,
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        magmoms::Ptr{Cdouble},
        tensor_rank::Cint,
        num_atom::Cint,
        with_time_reversal::Cint,
        is_axial::Cint,
        symprec::Cdouble,
    )::Cint
    check_error()
    rotations = map(
        SMatrix{3,3,Int32,9} ∘ transpose, eachslice(rotations[:, :, begin:num_sym]; dims=3)
    )  # Remember to transpose, see https://github.com/singularitti/Spglib.jl/blob/8aed6e0/src/core.jl#L195-L198
    translations = map(SVector{3,Float64}, eachcol(translations[:, begin:num_sym]))
    return rotations, translations, spin_flips[begin:num_sym]
end

struct SpglibMagneticDataset <: AbstractDataset
    uni_number::Cint
    msg_type::Cint
    hall_number::Cint
    tensor_rank::Cint
    n_operations::Cint
    rotations::Ptr{NTuple{3,NTuple{3,Cint}}}
    translations::Ptr{NTuple{3,Cdouble}}
    time_reversals::Ptr{Cint}
    n_atoms::Cint
    equivalent_atoms::Ptr{Cint}
    transformation_matrix::NTuple{3,NTuple{3,Cdouble}}
    origin_shift::NTuple{3,Cdouble}
    n_std_atoms::Cint
    std_lattice::NTuple{3,NTuple{3,Cdouble}}
    std_types::Ptr{Cint}
    std_positions::Ptr{NTuple{3,Cdouble}}
    std_tensors::Ptr{Cdouble}
    std_rotation_matrix::NTuple{3,NTuple{3,Cdouble}}
    primitive_lattice::NTuple{3,NTuple{3,Cdouble}}
end

"""
    MagneticDataset(uni_number, msg_type, hall_number, tensor_rank, n_operations, rotations, translations, time_reversals, n_atoms, equivalent_atoms, transformation_matrix, origin_shift, n_std_atoms, std_lattice, std_types, std_positions, std_tensors, std_rotation_matrix, primitive_lattice)

Represent `MagneticDataset`, see its [official documentation](https://spglib.readthedocs.io/en/latest/magnetic_dataset.html).

# Arguments
- `uni_number`: UNI number, from ``1`` to ``1651``.
- `msg_type`: Magnetic space groups (MSG) are classified by their family space group (FSG) and maximal space subgroup (XSG).

  The FSG is the nonmagnetic space group obtained by ignoring time-reversal in the MSG.
  The XSG is the space group obtained by removing time-reversal operations from the MSG.

  - msg_type == 1 (type I): MSG, XSG, and FSG are all isomorphic.
  - msg_type == 2 (type II): XSG and FSG are isomorphic; the MSG is generated from the XSG and pure time-reversal operations.
  - msg_type == 3 (type III): The XSG is a proper subgroup of the MSG with isomorphic translational subgroups.
  - msg_type == 4 (type IV): The XSG is a proper subgroup of the MSG with an isomorphic point group.
- `hall_number`: Hall number of the FSG (types I–III) or of the XSG (type IV).
- `tensor_rank`: Rank of magmoms.
- `n_operations`: Number of magnetic symmetry operations.
- `rotations`: Rotation matrices of the magnetic symmetry operations.
- `translations`: Translation vectors of the magnetic symmetry operations.
- `time_reversals`: Time-reversal flags for the magnetic symmetry operations. `true` indicates a time-reversal operation, and `false` indicates an ordinary operation.
- `n_atoms`: Number of atoms in the input cell.
- `equivalent_atoms`: Symmetrically equivalent atoms, where 'symmetrically equivalent' refers to the found symmetry operations.
- `transformation_matrix`: Transformation matrix (``3 × 3``) from the input lattice to the standardized lattice.
- `origin_shift`: Origin shift from the standardized origin to the input origin.
- `n_std_atoms`: Number of atoms in the standardized unit cell.
- `std_lattice`: Lattice vectors of the standardized unit cell.
- `std_types`: Identity numbers of atoms in the standardized unit cell.
- `std_positions`: Fractional coordinates of atoms in the standardized unit cell.
- `std_tensors`: Magnetic moments in the standardized unit cell: a length-`n_std_atoms` vector of scalars for collinear moments, or a length-`n_std_atoms` vector of 3-vectors for noncollinear moments.
- `std_rotation_matrix`: Rigid rotation matrix from the standardized basis vectors to an idealized, standardized, orthonormal basis.
- `primitive_lattice`: Basis vectors of a primitive cell.

See also [`get_magnetic_dataset`](@ref).
"""
@struct_hash_equal_isequal struct MagneticDataset <: AbstractDataset
    "UNI number, from ``1`` to ``1651``"
    uni_number::Int32
    """Magnetic space groups (MSG) are classified by their family space group (FSG) and
    maximal space subgroup (XSG).

    The FSG is the nonmagnetic space group obtained by ignoring time-reversal in the MSG.
    The XSG is the space group obtained by removing time-reversal operations from the MSG.

    - msg_type == 1 (type I):
        MSG, XSG, and FSG are all isomorphic.
    - msg_type == 2 (type II):
        XSG and FSG are isomorphic; the MSG is generated from the XSG and pure time-reversal operations.
    - msg_type == 3 (type III):
        The XSG is a proper subgroup of the MSG with isomorphic translational subgroups.
    - msg_type == 4 (type IV):
        The XSG is a proper subgroup of the MSG with an isomorphic point group.
    """
    msg_type::Int32
    "Hall number of the FSG (types I–III) or of the XSG (type IV)"
    hall_number::Int32
    "Rank of magmoms"
    tensor_rank::Int32
    "Number of magnetic symmetry operations"
    n_operations::Int32
    "Rotation matrices of the magnetic symmetry operations"
    rotations::Vector{SMatrix{3,3,Int32,9}}
    "Translation vectors of the magnetic symmetry operations"
    translations::Vector{SVector{3,Float64}}
    """Time-reversal flags for the magnetic symmetry operations.

    `true` indicates a time-reversal operation, and `false` indicates an ordinary operation.
    """
    time_reversals::BitVector
    "Number of atoms in the input cell"
    n_atoms::Int32
    "Symmetrically equivalent atoms, where 'symmetrically equivalent' refers to the found symmetry operations."
    equivalent_atoms::Vector{Int32}
    "Transformation matrix (``3 \times 3``) from the input lattice to the standardized lattice."
    transformation_matrix::SMatrix{3,3,Float64,9}
    "Origin shift from the standardized origin to the input origin."
    origin_shift::SVector{3,Float64}
    "Number of atoms in the standardized unit cell"
    n_std_atoms::Int32
    "Lattice vectors of the standardized unit cell"
    std_lattice::Lattice{Float64}
    "Identity numbers of atoms in the standardized unit cell"
    std_types::Vector{Int32}
    "Fractional coordinates of atoms in the standardized unit cell"
    std_positions::Vector{SVector{3,Float64}}
    """Magnetic moments in the standardized unit cell: a length-`n_std_atoms` vector
    of scalars for collinear moments, or a length-`n_std_atoms` vector of 3-vectors
    for noncollinear moments.
    """
    std_tensors::Vector{Union{Float64,SVector{3,Float64}}}
    "Rigid rotation matrix from the standardized basis vectors to an idealized, standardized, orthonormal basis."
    std_rotation_matrix::SMatrix{3,3,Float64,9}
    "Basis vectors of a primitive cell"
    primitive_lattice::Lattice{Float64}
end
function MagneticDataset(dataset::SpglibMagneticDataset)
    rotations =
        transpose.(
            _convert(SMatrix{3,3,Int32}, unsafe_load(dataset.rotations, i)) for
            i in Base.OneTo(dataset.n_operations)
        )
    translations =
        SVector{
            3
        }.(unsafe_load(dataset.translations, i) for i in Base.OneTo(dataset.n_operations))
    time_reversals =
        Bool.(unsafe_wrap(Vector{Int32}, dataset.time_reversals, dataset.n_operations))
    equivalent_atoms =  # Need to add 1 because of C-index starts from 0
        unsafe_wrap(Vector{Int32}, dataset.equivalent_atoms, dataset.n_atoms) .+ 1
    transformation_matrix = transpose(
        _convert(SMatrix{3,3,Float64}, dataset.transformation_matrix)
    )
    std_lattice = Lattice(transpose(_convert(SMatrix{3,3,Float64}, dataset.std_lattice)))
    std_types = unsafe_wrap(Vector{Int32}, dataset.std_types, dataset.n_std_atoms)
    std_positions =
        SVector{
            3
        }.(unsafe_load(dataset.std_positions, i) for i in Base.OneTo(dataset.n_std_atoms))
    std_tensors = if iszero(dataset.tensor_rank)  # Collinear spin
        unsafe_wrap(Vector{Float64}, dataset.std_tensors, dataset.n_std_atoms)
    else  # Non-collinear spin
        SVector{
            3
        }.(
            eachcol(
                unsafe_wrap(
                    Matrix{Float64},
                    dataset.std_tensors,
                    (3, Int(dataset.n_std_atoms)),  # Issue to Julia community
                ),
            ),
        )
    end
    std_rotation_matrix = transpose(
        _convert(SMatrix{3,3,Float64}, dataset.std_rotation_matrix)
    )
    primitive_lattice = Lattice(
        transpose(_convert(SMatrix{3,3,Float64}, dataset.primitive_lattice))
    )
    return MagneticDataset(
        dataset.uni_number,
        dataset.msg_type,
        dataset.hall_number,
        dataset.tensor_rank,
        dataset.n_operations,
        rotations,
        translations,
        time_reversals,
        dataset.n_atoms,
        equivalent_atoms,
        transformation_matrix,
        dataset.origin_shift,
        dataset.n_std_atoms,
        std_lattice,
        std_types,
        std_positions,
        std_tensors,
        std_rotation_matrix,
        primitive_lattice,
    )
end

"""
    get_magnetic_dataset(cell::SpglibCell, symprec=1e-5)

Return magnetic symmetry operations and standardized structure of given structure with site tensors.

The description of returned dataset is given at [Magnetic dataset (experimental)](@ref).
"""
function get_magnetic_dataset(cell::SpglibCell, symprec=1e-5)
    lattice, positions, atoms, magmoms = _unwrap_convert(cell)
    tensor_rank = ndims(magmoms) - 1  # See https://github.com/spglib/spglib/blob/v2.1.0/python/spglib/spglib.py#L275-L276 & https://github.com/spglib/spglib/blob/v2.1.0/python/spglib/spglib.py#L615
    is_axial = iszero(tensor_rank) ? false : true  # Collinear spin & non-collinear spin
    ptr = @ccall libsymspg.spg_get_magnetic_dataset(
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        magmoms::Ptr{Cdouble},
        tensor_rank::Cint,
        natoms(cell)::Cint,
        is_axial::Cint,
        symprec::Cdouble,
    )::Ptr{SpglibMagneticDataset}
    if ptr == C_NULL
        check_error()
    else
        raw_dataset = unsafe_load(ptr)
        return MagneticDataset(raw_dataset)
    end
end

# Python version: https://github.com/spglib/spglib/blob/v2.1.0/python/spglib/spglib.py#L1108-L1138
"""
    get_magnetic_symmetry_from_database(uni_number, hall_number=0)

Accesses magnetic space-group operations in the built-in database using UNI number (from ``1`` to ``1651``).

Optionally alternative settings can be specified with hall_number. For type-I, type-II, and
type-III magnetic space groups, `hall_number` changes settings in family space group. For
type-IV, `hall_number` changes settings in maximal space group. When `hall_number=0`, the
smallest hall number corresponding to `uni_number` is used.
"""
function get_magnetic_symmetry_from_database(uni_number, hall_number=0)
    @assert 1 <= uni_number <= 1651  # See https://github.com/spglib/spglib/blob/77a8e5d/src/spglib.h#L390
    _rotations = Array{Cint,3}(undef, 3, 3, 384)  # See https://github.com/spglib/spglib/blob/v2.1.0/src/spglib.h#L398
    _translations = Matrix{Cdouble}(undef, 3, 384)
    _time_reversals = Vector{Cint}(undef, 384)  # It cannot be `Vector{Bool}`!
    num_sym = @ccall libsymspg.spg_get_magnetic_symmetry_from_database(
        _rotations::Ptr{Cint},
        _translations::Ptr{Cdouble},
        _time_reversals::Ptr{Cint},
        uni_number::Cint,
        hall_number::Cint,
    )::Cint
    check_error()
    rotations = map(
        SMatrix{3,3,Int32,9} ∘ transpose, eachslice(_rotations[:, :, begin:num_sym]; dims=3)
    )  # Remember to transpose, see https://github.com/singularitti/Spglib.jl/blob/8aed6e0/src/core.jl#L195-L198
    translations = map(SVector{3,Float64}, eachcol(_translations[:, begin:num_sym]))
    time_reversals = Bool.(_time_reversals[begin:num_sym])
    return rotations, translations, time_reversals
end

struct SpglibMagneticSpacegroupType <: AbstractSpacegroupType
    uni_number::Cint
    litvin_number::Cint
    bns_number::NTuple{8,Cchar}
    og_number::NTuple{12,Cchar}
    number::Cint
    type::Cint
end

"""
   MagneticSpacegroupType(uni_number, litvin_number, bns_number, og_number, number, type)

Represent `SpglibMagneticSpacegroupType`, see its [official documentation](https://spglib.readthedocs.io/en/latest/api.html#spg-get-magnetic-spacegroup-type).

# Arguments
- `uni_number::Int32`: Serial number of UNI (or BNS) symbols.
- `litvin_number::Int32`: Serial number in Litvin's [Magnetic Group Tables](https://www.iucr.org/publ/978-0-9553602-2-0).
- `bns_number::String`: BNS number, e.g. `"151.32"`.
- `og_number::String`: OG number, e.g. `"153.4.1270"`.
- `number::Int32`: ITA's serial number of space group for reference setting.
- `type::Int32`: Type of MSG, from ``1`` to ``4``.

See also [`get_magnetic_spacegroup_type`](@ref), [`get_magnetic_spacegroup_type_from_symmetry`](@ref).
"""
struct MagneticSpacegroupType <: AbstractSpacegroupType
    uni_number::Int32
    litvin_number::Int32
    bns_number::String
    og_number::String
    number::Int32
    type::Int32
end
function MagneticSpacegroupType(spgtype::SpglibMagneticSpacegroupType)
    bns_number = tostring(spgtype.bns_number)
    og_number = tostring(spgtype.og_number)
    return MagneticSpacegroupType(
        spgtype.uni_number,
        spgtype.litvin_number,
        bns_number,
        og_number,
        spgtype.number,
        spgtype.type,
    )
end

"""
    get_magnetic_spacegroup_type(uni_number)

Accesses magnetic space-group type by serial number (from ``1`` to ``1651``) of UNI (or BNS) symbols.
"""
function get_magnetic_spacegroup_type(uni_number)
    spgtype = @ccall libsymspg.spg_get_magnetic_spacegroup_type(
        uni_number::Cint
    )::SpglibMagneticSpacegroupType
    check_error()
    return MagneticSpacegroupType(spgtype)
end

"""
    get_magnetic_spacegroup_type_from_symmetry(rotations, translations, time_reversals, lattice::Lattice, symprec=1e-5)

Determine magnetic space-group type from magnetic symmetry operations.

`time_reversals` takes `false` for ordinary operations and `true` for time-reversal operations.
"""
function get_magnetic_spacegroup_type_from_symmetry(
    rotations, translations, time_reversals, lattice::Lattice, symprec=1e-5
)
    if length(rotations) != length(translations)
        throw(DimensionMismatch("the numbers of rotations and translations are different!"))
    end
    num_sym = length(translations)
    rotations = Base.cconvert(Array{Cint,3}, cat(transpose.(rotations)...; dims=3))
    translations = Base.cconvert(Matrix{Cdouble}, reduce(hcat, translations))
    time_reversals = Base.cconvert(Vector{Cint}, time_reversals)
    lattice = Base.cconvert(Matrix{Cdouble}, transpose(lattice))   # `transpose` must before `cconvert`!
    spgtype = @ccall libsymspg.spg_get_magnetic_spacegroup_type_from_symmetry(
        rotations::Ptr{Cint},
        translations::Ptr{Cdouble},
        time_reversals::Ptr{Cint},
        num_sym::Cint,
        lattice::Ptr{Cdouble},
        symprec::Cdouble,
    )::SpglibMagneticSpacegroupType
    check_error()
    return MagneticSpacegroupType(spgtype)
end
