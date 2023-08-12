export MagneticDataset

function get_symmetry_with_collinear_spin(cell::SpglibCell, symprec=1e-5)
    lattice, positions, atoms, magmoms = _expand_cell(cell)
    n = length(cell.magmoms)
    # See https://github.com/spglib/spglib/blob/42527b0/python/spglib/spglib.py#L270
    max_size = 96n  # 96 = 48 × 2 since we have spins
    rotations = Array{Cint,3}(undef, 3, 3, max_size)
    translations = Matrix{Cdouble}(undef, 3, max_size)
    equivalent_atoms = Vector{Cint}(undef, n)
    nsym = @ccall libsymspg.spg_get_symmetry_with_collinear_spin(
        rotations::Ptr{Cint},
        translations::Ptr{Cdouble},
        equivalent_atoms::Ptr{Cint},
        max_size::Cint,
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        magmoms::Ptr{Cdouble},
        n::Cint,
        symprec::Cdouble,
    )::Cint
    check_error()
    rotations, translations = map(
        SMatrix{3,3,Int32,9}, eachslice(rotations[:, :, 1:nsym]; dims=3)
    ),
    map(SVector{3,Float64}, eachcol(translations[:, 1:nsym]))
    return rotations, translations
end
const get_magnetic_symmetry = get_symmetry_with_collinear_spin

function get_symmetry_with_site_tensors(
    cell::SpglibCell, symprec=1e-5; with_time_reversal=true, is_axial=false
)
    lattice, positions, atoms, magmoms = _expand_cell(cell)
    n = length(cell.magmoms)
    # See https://github.com/spglib/spglib/blob/42527b0/python/spglib/spglib.py#L270
    max_size = 96n  # 96 = 48 × 2 since we have spins
    rotations = Array{Cint,3}(undef, 3, 3, max_size)
    translations = Matrix{Cdouble}(undef, 3, max_size)
    equivalent_atoms = Vector{Cint}(undef, n)
    primitive_lattice = zeros(Cdouble, 3, 3)
    spin_flips = if ndims(magmoms) == 1
        zeros(length(rotations))
    else
        nothing
    end
    nsym = @ccall libsymspg.spg_get_symmetry_with_site_tensors(
        rotations::Ptr{Cint},
        translations::Ptr{Cdouble},
        equivalent_atoms::Ptr{Cint},
        primitive_lattice::Ptr{Cdouble},
        spin_flips::Ptr{Cint},
        max_size::Cint,
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        tensors::Ptr{Cdouble},
        tensor_rank::Cint,
        natoms(cell)::Cint,
        with_time_reversal::Cint,
        is_axial::Cint,
        symprec::Cdouble,
    )::Cint
    check_error()
    rotations, translations = map(
        SMatrix{3,3,Int32,9}, eachslice(rotations[:, :, 1:nsym]; dims=3)
    ),
    map(SVector{3,Float64}, eachcol(translations[:, 1:nsym]))
    return rotations, translations
end

struct SpglibMagneticDataset
    uni_number::Cint
    msg_type::Cint
    hall_number::Cint
    tensor_rank::Cint
    n_operations::Cint
    rotations::Ptr{Cint}
    translations::Ptr{Cdouble}
    time_reversals::Ptr{Cint}
    n_atoms::Cint
    equivalent_atoms::Ptr{Cint}
    transformation_matrix::NTuple{9,Cdouble}
    origin_shift::NTuple{3,Cdouble}
    n_std_atoms::Cint
    std_lattice::NTuple{9,Cdouble}
    std_types::Ptr{Cint}
    std_positions::Ptr{Cdouble}
    std_tensors::Ptr{Cdouble}
    std_rotation_matrix::NTuple{9,Cdouble}
    primitive_lattice::NTuple{9,Cdouble}
end

struct MagneticDataset
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
    std_tensors::Vector{Float64}
    std_rotation_matrix::SMatrix{3,3,Float64,9}
    primitive_lattice::Lattice{Float64}
end

function get_magnetic_dataset(
    cell::SpglibCell, tensor_rank::Cint, is_axial=false, symprec=1e-5
)
    lattice, positions, atoms, magmoms = _expand_cell(cell)
    ptr = @ccall libsymspg.spg_get_magnetic_dataset(
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        tensors::Ptr{Cdouble},
        tensor_rank::Cint,
        natoms(cell)::Cint,
        is_axial::Cint,
        symprec::Cdouble,
    )::Ptr{SpglibMagneticDataset}
    if ptr == C_NULL
        check_error()
    else
        raw = unsafe_load(ptr)
        return convert(Dataset, raw)
    end
end

function get_magnetic_symmetry_from_database(
    cell::SpglibCell, uni_number::Cint, hall_number::Cint
)
    @assert 1 <= uni_number <= 1651  # See https://github.com/spglib/spglib/blob/77a8e5d/src/spglib.h#L390
    @ccall libsymspg.spg_get_magnetic_symmetry_from_database(
        rotations::Ptr{Cint},
        translations::Ptr{Cdouble},
        time_reversals::Ptr{Cint},
        uni_number::Cint,
        hall_number::Cint,
    )::Cint
end

struct SpglibMagneticSpacegroupType
    uni_number::Cint
    litvin_number::Cint
    bns_number::NTuple{8,Cchar}
    og_number::NTuple{12,Cchar}
    number::Cint
    type::Cint
end

struct MagneticSpacegroupType
    uni_number::Int32
    litvin_number::Int32
    bns_number::String
    og_number::String
    number::Int32
    type::Int32
end

function get_magnetic_spacegroup_type(uni_number::Integer)
    spgtype = @ccall libsymspg.spg_get_magnetic_spacegroup_type(
        uni_number::Cint
    )::SpglibMagneticSpacegroupType
    check_error()
    return convert(MagneticSpacegroupType, spgtype)
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
    return convert(MagneticSpacegroupType, spgtype)
end

function Base.convert(::Type{MagneticDataset}, dataset::SpglibMagneticDataset)
    rotations = [
        _convert(SMatrix{3,3,Int32}, unsafe_load(dataset.rotations, i)) for
        i in Base.OneTo(dataset.n_operations)
    ]
    translations = [
        SVector{3}(unsafe_load(dataset.translations, i)) for
        i in Base.OneTo(dataset.n_operations)
    ]
    time_reversals = unsafe_wrap(
        Vector{Int32}, dataset.time_reversals, dataset.n_operations
    )
    equivalent_atoms = unsafe_wrap(Vector{Int32}, dataset.equivalent_atoms, dataset.n_atoms)
    std_lattice = Lattice(transpose(_convert(SMatrix{3,3,Float64}, dataset.std_lattice)))
    std_types = unsafe_wrap(Vector{Int32}, dataset.std_types, dataset.n_std_atoms)
    std_positions = [
        SVector{3}(unsafe_load(dataset.std_positions, i)) for
        i in Base.OneTo(dataset.n_std_atoms)
    ]
    std_tensors = unsafe_wrap(Vector{Float64}, dataset.std_tensors, dataset.n_std_atoms)
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
        dataset.transformation_matrix,
        dataset.origin_shift,
        dataset.n_std_atoms,
        std_lattice,
        std_types,
        std_positions,
        std_tensors,
        dataset.std_rotation_matrix,
        dataset.primitive_lattice,
    )
end
