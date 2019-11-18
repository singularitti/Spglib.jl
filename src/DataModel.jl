"""
# module DataModel



# Examples

```jldoctest
julia>
```
"""
module DataModel

using Parameters

export Cell, Dataset, SpacegroupType

@with_kw struct Cell{L<:AbstractMatrix,P<:AbstractMatrix,N<:AbstractVector,M<:Union{AbstractVector,Nothing}}
    lattice::L
    positions::P
    numbers::N
    magmoms::M = nothing
end
Cell(lattice, positions, numbers) = Cell(lattice, positions, numbers, nothing)

struct Dataset
    spacegroup_number::Cint
    hall_number::Cint
    international_symbol::NTuple{11,Cchar}
    hall_symbol::NTuple{17,Cchar}
    choice::NTuple{6,Cchar}
    transformation_matrix::Matrix{Cdouble}
    origin_shift::Vector{Cdouble}
    n_operations::Cint
    rotations::Ptr{Cint}
    translations::Ptr{Cdouble}
    n_atoms::Cint
    wyckoffs::Ptr{Cint}
    site_symmetry_symbols::Ptr{Cchar}
    equivalent_atoms::Ptr{Cint}
    mapping_to_primitive::Ptr{Cint}
    n_std_atoms::Cint
    std_lattice::Matrix{Cdouble}
    std_types::Ptr{Cint}
    std_positions::Ptr{Cdouble}
    std_rotation_matrix::Matrix{Cdouble}
    std_mapping_to_primitive::Ptr{Cint}
    pointgroup_symbol::NTuple{6,Cchar}
end

struct SpacegroupType
    number::Cint
    international_short::Cstring
    international_full::Cstring
    international::Cstring
    schoenflies::Cstring
    hall_symbol::Cstring
    choice::Cstring
    pointgroup_international::Cstring
    pointgroup_schoenflies::Cstring
    arithmetic_crystal_class_number::Cint
    arithmetic_crystal_class_symbol::Cstring
end

end
