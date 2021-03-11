module DataModel

export Cell, Dataset, SpaceGroup

struct Cell{
    L<:AbstractVecOrMat,
    P<:AbstractVecOrMat,
    N<:AbstractVector,
    M<:Union{AbstractVector,Nothing},
}
    lattice::L
    positions::P
    numbers::N
    magmoms::M
end
Cell(lattice, positions, numbers) = Cell(lattice, positions, numbers, nothing)

# This is an internal type, do not export!
struct Cdataset
    spacegroup_number::Cint
    hall_number::Cint
    international_symbol::NTuple{11,UInt8}
    hall_symbol::NTuple{17,UInt8}
    choice::NTuple{6,UInt8}
    transformation_matrix::NTuple{3,NTuple{3,Cdouble}}
    origin_shift::NTuple{3,Cdouble}
    n_operations::Cint
    rotations::Ptr{Cvoid}
    translations::Ptr{Cvoid}
    n_atoms::Cint
    wyckoffs::Ptr{Cint}
    equivalent_atoms::Ptr{Cint}
    n_std_atoms::Cint
    std_lattice::NTuple{3,NTuple{3,Cdouble}}
    std_types::Ptr{Cint}
    std_positions::Ptr{Cvoid}
    pointgroup_symbol::NTuple{6,UInt8}
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
    rotations::Vector{Matrix{Float64}}
    translations::Vector{Float64}
    n_atoms::Int
    wyckoffs::Vector{Int}
    equivalent_atoms::Vector{Int}
    n_std_atoms::Int
    std_lattice::Matrix{Float64}
    std_types::Vector{Int}
    std_positions::Matrix{Float64}
    pointgroup_symbol::String
end

# This is an internal type, do not export!
struct CspaceGroup
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

end
