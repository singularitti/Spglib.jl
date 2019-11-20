module SpgLib

using Parameters: @with_kw

export Cell, Dataset, SpaceGroup

@with_kw struct Cell{L<:AbstractVecOrMat,P<:AbstractVecOrMat,N<:AbstractVector,M<:Union{AbstractVector,Nothing}}
    lattice::L
    positions::P
    numbers::N
    magmoms::M = nothing
end
Cell(lattice, positions, numbers) = Cell(lattice, positions, numbers, nothing)

struct Dataset
    spacegroup_number::Int
    hall_number::Int
    international_symbol::String
    hall_symbol::String
    choice::String
    transformation_matrix::Matrix{Float64}
    origin_shift::Vector{Float64}
    n_operations::Int
    rotations::Ptr{Cvoid}
    translations::Ptr{Cvoid}
    n_atoms::Int
    wyckoffs::Ptr{Int}
    equivalent_atoms::Ptr{Int}
    n_std_atoms::Int
    std_lattice::Matrix{Float64}
    std_types::Ptr{Int}
    std_positions::Ptr{Cvoid}
    pointgroup_symbol::String
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

include("Wrapper/Wrapper.jl")
include("FFI.jl")

end # module
