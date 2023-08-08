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
