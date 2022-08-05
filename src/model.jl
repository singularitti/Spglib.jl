using StaticArrays: MMatrix, MVector
using StructHelpers: @batteries

export Cell, Dataset, SpacegroupType, basis_vectors

"""
    Cell(lattice, positions, types, magmoms=zeros(length(types)))

The basic input data type of `Spglib`.

Lattice parameters `lattice` are given by a ``3×3`` matrix with floating point values,
where ``𝐚``, ``𝐛``, and ``𝐜`` are given as columns.
Fractional atomic positions `positions` are given
by a ``3×N`` matrix with floating point values, where ``N`` is the number of atoms.
Numbers to distinguish atomic species `types` are given by a list of ``N`` integers.
The collinear polarizations `magmoms` only work with `get_symmetry` and are given
as a list of ``N`` floating point values.
"""
struct Cell{L,P,T,M}
    lattice::MMatrix{3,3,L,9}
    positions::Vector{MVector{3,P}}
    types::Vector{T}
    magmoms::Vector{M}
end
function Cell(lattice, positions, types, magmoms = zeros(length(types)))
    if lattice isa AbstractVector
        lattice = hcat(lattice...)
    end
    N = length(types)
    if positions isa AbstractMatrix
        P = eltype(positions)
        if size(positions) == (3, N)
            positions = collect(eachcol(positions))
        elseif size(positions) == (N, 3)
            positions = collect(eachrow(positions))
        else
            #! format: off
            throw(DimensionMismatch( "the `positions` has a different number of atoms from the `types`!"))
            #! format: on
        end
    else  # positions isa AbstractVector or a Tuple
        P = eltype(Base.promote_typeof(positions...))
        positions = collect(map(MVector{3,P}, positions))
    end
    L, T, M = eltype(lattice), eltype(types), eltype(magmoms)
    return Cell{L,P,T,M}(lattice, positions, types, magmoms)
end

@batteries Cell eq = true hash = true

"""
    basis_vectors(cell::Cell)

Return the three basis vectors from `cell`.
"""
function basis_vectors(cell::Cell)
    lattice = cell.lattice
    return lattice[:, 1], lattice[:, 2], lattice[:, 3]
end

# This is an internal function, do not export!
function _expand_cell(cell::Cell)
    @unpack lattice, positions, types, magmoms = cell
    # Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L32-L35 and https://github.com/spglib/spglib/blob/444e061/python/spglib/spglib.py#L953-L975
    clattice = Base.cconvert(Matrix{Cdouble}, lattice) |> transpose
    cpositions = Base.cconvert(Matrix{Cdouble}, positions)
    ctypes = Cint[findfirst(isequal(u), unique(types)) for u in types]
    if magmoms !== nothing
        magmoms = Base.cconvert(Vector{Cdouble}, magmoms)
    end
    return Cell(clattice, cpositions, ctypes, magmoms)
end

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

"""
    SpglibSpacegroupType(number, international_short, international_full, international, schoenflies, hall_symbol, choice, pointgroup_international, pointgroup_schoenflies, arithmetic_crystal_class_number, arithmetic_crystal_class_symbol)

Represent `SpglibSpacegroupType`, see its [official documentation](https://spglib.github.io/spglib/api.html#spg-get-spacegroup-type).
"""
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

struct SpglibError <: Exception
    msg::AbstractString
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

@static if get_version() >= v"1.15.0"
    include("new.jl")
else
    include("old.jl")
end

function Base.show(io::IO, cell::Cell{N}) where {N}
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(cell)
        Base.show_default(IOContext(io, :limit => true), cell)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, string(typeof(cell)))
        println(io, "lattice:")
        print(io, " ")
        show(io, "text/plain", cell.lattice)
        println(io, "")
        println(io, "$N atomic positions:")
        print(io, " ")
        show(io, "text/plain", cell.positions)
        println(io, "")
        println(io, "$N atoms:")
        println(io, " ", cell.types)
    end
end
