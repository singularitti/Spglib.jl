export Cell, Dataset, SpacegroupType

"""
    Cell(lattice, positions, types, magmoms=nothing)

The basic input data type of `Spglib`.

Lattice parameters `lattice` are given by a 3×3 matrix with floating point values,
where 𝐚, 𝐛, and 𝐜 are given as rows, which results in the transpose of the
definition for C-API. Fractional atomic positions `positions` are given
by a N×3 matrix with floating point values, where N is the number of atoms.
Numbers to distinguish atomic species `types` are given by a list of N integers.
The collinear polarizations `magmoms` only work with `get_symmetry` and are given
as a list of N floating point values.
"""
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
struct SpglibSpacegroupType
    number::Cint
    international_short::NTuple{11,Cchar}
    international_full::NTuple{20,Cchar}
    international::NTuple{32,Cchar}
    schoenflies::NTuple{7,Cchar}
    hall_symbol::NTuple{17,Cchar}
    choice::NTuple{6,Cchar}
    pointgroup_international::NTuple{6,Cchar}
    pointgroup_schoenflies::NTuple{4,Cchar}
    arithmetic_crystal_class_number::Cint
    arithmetic_crystal_class_symbol::NTuple{7,Cchar}
end

"This represents `SpglibSpacegroupType`, see https://spglib.github.io/spglib/api.html#spg-get-spacegroup-type."
struct SpacegroupType
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

"""
    get_version()

Obtain the version number of `spglib`.

This is the mergence of `spg_get_major_version`, `spg_get_minor_version`, and `spg_get_micro_version` in its C-API.
"""
function get_version()
    major = ccall((:spg_get_major_version, libsymspg), Cint, ())
    minor = ccall((:spg_get_minor_version, libsymspg), Cint, ())
    micro = ccall((:spg_get_micro_version, libsymspg), Cint, ())
    return VersionNumber(major, minor, micro)
end

tuple2matrix(t::NTuple{9}) = hcat(Iterators.partition(t, 3)...)

function rotsFromTuple(rotsTuple::AbstractVector{NTuple{9,Int32}}, nop::Integer)
    r = Array{Int64,3}(undef, 3, 3, nop)
    for i in 1:nop
        r[:, :, i] = reshape(collect(rotsTuple[i]), 3, 3)
    end
    return r
end

function transFromTuple(transTuple::AbstractVector{NTuple{3,Float64}}, nop::Integer)
    t = Matrix{Float64}(undef, 3, nop)
    for i in 1:nop
        t[:, i] = collect(transTuple[i])
    end
    return t
end

const LETTERS = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

if get_version() >= v"1.15.0"
    include("new.jl")
else
    include("old.jl")
end
