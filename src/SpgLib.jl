module SpgLib

using Parameters: @with_kw

export Cell, SpaceGroup

@with_kw struct Cell{L<:AbstractVecOrMat,P<:AbstractVecOrMat,N<:AbstractVector,M<:Union{AbstractVector,Nothing}}
    lattice::L
    positions::P
    numbers::N
    magmoms::M = nothing
end
Cell(lattice, positions, numbers) = Cell(lattice, positions, numbers, nothing)

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
