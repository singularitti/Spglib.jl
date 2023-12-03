export MagneticDataset,
    MagneticSpacegroupType,
    get_symmetry_with_collinear_spin,
    get_symmetry_with_site_tensors,
    get_magnetic_symmetry,
    get_magnetic_dataset,
    get_magnetic_symmetry_from_database,
    get_magnetic_spacegroup_type

# Python version: https://github.com/spglib/spglib/blob/42527b0/python/spglib/spglib.py#L182-L319
function get_symmetry_with_collinear_spin(cell::SpglibCell, symprec=1e-5)
    lattice, positions, atoms, spins = _unwrap_convert(cell)
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
        spins::Ptr{Cdouble},
        num_atom::Cint,
        symprec::Cdouble,
    )::Cint
    check_error()
    rotations = map(
        SMatrix{3,3,Int32,9} ∘ transpose, eachslice(rotations[:, :, 1:num_sym]; dims=3)
    )  # Remember to transpose, see https://github.com/singularitti/Spglib.jl/blob/8aed6e0/src/core.jl#L195-L198
    translations = map(SVector{3,Float64}, eachcol(translations[:, 1:num_sym]))
    return rotations, translations, equivalent_atoms .+ 1
end
const get_magnetic_symmetry = get_symmetry_with_collinear_spin

function get_symmetry_with_site_tensors(
    cell::SpglibCell, symprec=1e-5; with_time_reversal=true, is_axial=false
)
    lattice, positions, atoms, magmoms = _unwrap_convert(cell)
    num_atom = natoms(cell)
    # See https://github.com/spglib/spglib/blob/42527b0/python/spglib/spglib.py#L270
    max_size = 96num_atom  # 96 = 48 × 2 since we have spins
    rotations = Array{Cint,3}(undef, 3, 3, max_size)
    translations = Matrix{Cdouble}(undef, 3, max_size)
    equivalent_atoms = Vector{Cint}(undef, num_atom)
    primitive_lattice = zeros(Cdouble, 3, 3)
    spin_flips = if isone(ndims(magmoms))
        zeros(Cint, length(rotations))
    else
        nothing
    end
    tensor_rank = ndims(magmoms) - 1  # See https://github.com/spglib/spglib/blob/v2.1.0/python/spglib/spglib.py#L275-L276 & https://github.com/spglib/spglib/blob/v2.1.0/python/spglib/spglib.py#L615
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
        SMatrix{3,3,Int32,9} ∘ transpose, eachslice(rotations[:, :, 1:num_sym]; dims=3)
    )  # Remember to transpose, see https://github.com/singularitti/Spglib.jl/blob/8aed6e0/src/core.jl#L195-L198
    translations = map(SVector{3,Float64}, eachcol(translations[:, 1:num_sym]))
    return rotations, translations, spin_flips[1:num_sym]
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

@struct_hash_equal_isequal struct MagneticDataset <: AbstractDataset
    uni_number::Int32
    msg_type::Int32
    hall_number::Int32
    tensor_rank::Int32
    n_operations::Int32
    rotations::Vector{SMatrix{3,3,Int32,9}}
    translations::Vector{SVector{3,Float64}}
    time_reversals::Vector{Int32}
    n_atoms::Int32
    equivalent_atoms::Vector{Int32}
    transformation_matrix::SMatrix{3,3,Float64,9}
    origin_shift::SVector{3,Float64}
    n_std_atoms::Int32
    std_lattice::Lattice{Float64}
    std_types::Vector{Int32}
    std_positions::Vector{SVector{3,Float64}}
    std_tensors::Vector{Union{Float64,SVector{3,Float64}}}
    std_rotation_matrix::SMatrix{3,3,Float64,9}
    primitive_lattice::Lattice{Float64}
end
function MagneticDataset(dataset::SpglibMagneticDataset)
    rotations = transpose.(
        _convert(SMatrix{3,3,Int32}, unsafe_load(dataset.rotations, i)) for
        i in Base.OneTo(dataset.n_operations)
    )
    translations = SVector{3}.(
        unsafe_load(dataset.translations, i) for i in Base.OneTo(dataset.n_operations)
    )
    time_reversals = unsafe_wrap(
        Vector{Int32}, dataset.time_reversals, dataset.n_operations
    )
    equivalent_atoms =  # Need to add 1 because of C-index starts from 0
        unsafe_wrap(Vector{Int32}, dataset.equivalent_atoms, dataset.n_atoms) .+ 1
    transformation_matrix = transpose(
        _convert(SMatrix{3,3,Float64}, dataset.transformation_matrix)
    )
    std_lattice = Lattice(transpose(_convert(SMatrix{3,3,Float64}, dataset.std_lattice)))
    std_types = unsafe_wrap(Vector{Int32}, dataset.std_types, dataset.n_std_atoms)
    std_positions = SVector{3}.(
        unsafe_load(dataset.std_positions, i) for i in Base.OneTo(dataset.n_std_atoms)
    )
    std_tensors = if iszero(dataset.tensor_rank)  # Collinear spin
        unsafe_wrap(Vector{Float64}, dataset.std_tensors, dataset.n_std_atoms)
    else  # Non-collinear spin
        SVector{3}.(
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

function get_magnetic_spacegroup_type(uni_number)
    spgtype = @ccall libsymspg.spg_get_magnetic_spacegroup_type(
        uni_number::Cint
    )::SpglibMagneticSpacegroupType
    check_error()
    return MagneticSpacegroupType(spgtype)
end

function get_magnetic_spacegroup_type_from_symmetry(cell::SpglibCell, symprec=1e-5)
    rotations, translations = get_symmetry(cell, symprec)
    nsym = length(translations)
    rotations, translations = reduce(hcat, rotations), reduce(hcat, translations)
    time_reversals = zeros(Int32, nsym)
    spgtype = @ccall libsymspg.spg_get_magnetic_spacegroup_type_from_symmetry(
        rotations::Ptr{Cint},
        translations::Ptr{Cdouble},
        time_reversals::Ptr{Cint},
        nsym::Cint,
        Lattice(cell)::Ptr{Cdouble},
        symprec::Cdouble,
    )::SpglibMagneticSpacegroupType
    check_error()
    return MagneticSpacegroupType(spgtype)
end
