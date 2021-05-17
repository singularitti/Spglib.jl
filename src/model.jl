export Cell, Dataset, SpaceGroup

struct Cell{
    L<:AbstractVecOrMat,
    P<:AbstractVecOrMat,
    N<:AbstractVector,
    M<:Union{AbstractVector,Nothing},
}
    lattice::L
    positions::P
    types::N
    magmoms::M
end
Cell(lattice, positions, types) = Cell(lattice, positions, types, nothing)

# This is an internal type, do not export!
struct SpglibDataset
    spacegroup_number::Cint
    hall_number::Cint
    international_symbol::NTuple{11,Cchar}
    hall_symbol::NTuple{17,Cchar}
    choice::NTuple{6,Cchar}
    transformation_matrix::NTuple{9,Cdouble}
    origin_shift::NTuple{3,Cdouble}
    n_operations::Cint
    rotations::Ptr{NTuple{9,Cint}}
    translations::Ptr{NTuple{3,Cdouble}}
    n_atoms::Cint
    wyckoffs::Ptr{Cint}
    site_symmetry_symbols::Ptr{Tuple{7,Cchar}}
    equivalent_atoms::Ptr{Cint}
    crystallographic_orbits::Ptr{Cint}  # Added in v1.15.0
    primitive_lattice::NTuple{9,Cdouble}  # Added in v1.15.0
    mapping_to_primitive::Ptr{Cint}
    n_std_atoms::Cint
    std_lattice::NTuple{9,Cdouble}
    std_types::Ptr{Cint}
    std_positions::Ptr{NTuple{3,Cdouble}}
    std_rotation_matrix::NTuple{9,Cdouble}
    std_mapping_to_primitive::Ptr{Cint}
    pointgroup_symbol::NTuple{6,Cchar}
end

struct Dataset
    spacegroup_number::Int
    hall_number::Int
    international_symbol::String
    hall_symbol::String
    choice::String
    transformation_matrix::Matrix{Float64}
    origin_shift::Vector{Float64}
    n_operations::Int
    rotations::Matrix{Float64}
    translations::Vector{Float64}
    n_atoms::Int
    wyckoffs::Int
    site_symmetry_symbols::String
    equivalent_atoms::Int
    crystallographic_orbits::Int
    primitive_lattice::Matrix{Float64}
    mapping_to_primitive::Int
    n_std_atoms::Int
    std_lattice::Matrix{Float64}
    std_types::Int
    std_positions::Vector{Float64}
    std_rotation_matrix::Matrix{Float64}
    std_mapping_to_primitive::Int
    pointgroup_symbol::String
end

# This is an internal type, do not export!
struct SpglibSpacegroupType
    number::Cint
    international_short::NTuple{11,UInt8}
    international_full::NTuple{20,UInt8}
    international::NTuple{32,UInt8}
    schoenflies::NTuple{7,UInt8}
    hall_symbol::NTuple{17,UInt8}
    choice::NTuple{6,UInt8}
    pointgroup_international::NTuple{6,UInt8}
    pointgroup_schoenflies::NTuple{4,UInt8}
    arithmetic_crystal_class_number::Cint
    arithmetic_crystal_class_symbol::NTuple{7,UInt8}
end

struct SpaceGroup
    number::Int
    international_short::String
    international_full::String
    international::String
    schoenflies::String
    hall_symbol::String
    choice::String
    pointgroup_international::String
    pointgroup_schoenflies::String
    arithmetic_crystal_class_number::Int
    arithmetic_crystal_class_symbol::String
end
