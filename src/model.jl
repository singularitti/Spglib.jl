using StaticArrays: MMatrix, MVector
using StructEquality: @struct_hash_equal

export Cell, Dataset, SpacegroupType, basis_vectors, natoms

"""
    Cell(lattice, positions, types, magmoms=zeros(length(types)))

The basic input data type of `Spglib`.

Lattice parameters `lattice` are given by a ``3×3`` matrix with floating point values,
where ``𝐚``, ``𝐛``, and ``𝐜`` are given as columns.
Fractional atomic positions `positions` are given
by a vector of ``N`` vectors with floating point values, where ``N`` is the number of atoms.
Numbers to distinguish atomic species `types` are given by a list of ``N`` integers.
The collinear polarizations `magmoms` only work with `get_symmetry` and are given
as a list of ``N`` floating point values, or a vector of vectors.
"""
@struct_hash_equal struct Cell{L,P,T,M}
    lattice::MMatrix{3,3,L,9}
    positions::Vector{MVector{3,P}}
    types::Vector{T}
    magmoms::M
end
function Cell(lattice, positions, types, magmoms=nothing)
    if !(lattice isa AbstractMatrix)
        lattice = reduce(hcat, lattice)  # Use `reduce` can make it type stable
    end
    N = length(types)
    if positions isa AbstractMatrix
        P = eltype(positions)
        if size(positions) == (3, 3)
            error("ambiguous `positions` size 3×3! Use a vector of `Vector`s instead!")
        elseif size(positions) == (3, N)
            positions = collect(eachcol(positions))
        elseif size(positions) == (N, 3)
            positions = collect(eachrow(positions))
        else
            throw(
                DimensionMismatch(
                    "the `positions` has a different number of atoms from the `types`!"
                ),
            )
        end
    else  # positions isa AbstractVector or a Tuple
        P = eltype(Base.promote_typeof(positions...))
        positions = collect(map(MVector{3,P}, positions))
    end
    L, T, M = eltype(lattice), eltype(types), typeof(magmoms)
    return Cell{L,P,T,M}(lattice, positions, types, magmoms)
end

natoms(cell::Cell) = length(cell.types)

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
    lattice, positions, types, magmoms = cell.lattice,
    cell.positions, cell.types,
    cell.magmoms
    # Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L32-L35 and https://github.com/spglib/spglib/blob/444e061/python/spglib/spglib.py#L953-L975
    clattice = Base.cconvert(Matrix{Cdouble}, transpose(lattice))
    cpositions = Base.cconvert(Matrix{Cdouble}, reduce(hcat, positions))
    ctypes = Cint[findfirst(isequal(u), unique(types)) for u in types]
    if magmoms !== nothing
        magmoms = Base.cconvert(Vector{Cdouble}, magmoms)
    end
    return clattice, cpositions, ctypes, magmoms
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
    site_symmetry_symbols::Ptr{NTuple{7,Cchar}}
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

"""
    Dataset(spacegroup_number, hall_number, international_symbol, hall_symbol, choice, transformation_matrix, origin_shift, n_operations, rotations, translations, n_atoms, wyckoffs, site_symmetry_symbols, equivalent_atoms, crystallographic_orbits, primitive_lattice, mapping_to_primitive, n_std_atoms, std_lattice, std_types, std_positions, std_rotation_matrix, std_mapping_to_primitive, pointgroup_symbol)

Represent `SpglibDataset`, see its [official documentation](https://spglib.github.io/spglib/dataset.html#spglib-dataset).

!!! note
    Fields `crystallographic_orbits` and `primitive_lattice` are added after `spglib` `v1.15.0`.
"""
struct Dataset
    spacegroup_number::Int
    hall_number::Int
    international_symbol::String
    hall_symbol::String
    choice::String
    transformation_matrix::Matrix{Float64}
    origin_shift::Vector{Float64}
    n_operations::Int
    rotations::Array{Float64,3}
    translations::Matrix{Float64}
    n_atoms::Int
    wyckoffs::Vector{Char}
    site_symmetry_symbols::Vector{String}
    equivalent_atoms::Vector{Int}
    crystallographic_orbits::Vector{Int}
    primitive_lattice::Matrix{Float64}
    mapping_to_primitive::Vector{Int}
    n_std_atoms::Int
    std_lattice::Matrix{Float64}
    std_types::Vector{Int}
    std_positions::Matrix{Float64}
    std_rotation_matrix::Matrix{Float64}
    std_mapping_to_primitive::Vector{Int}
    pointgroup_symbol::String
end

function Base.convert(::Type{Dataset}, dataset::SpglibDataset)
    r = unsafe_wrap(Vector{NTuple{9,Cint}}, dataset.rotations, dataset.n_operations)
    t = unsafe_wrap(Vector{NTuple{3,Float64}}, dataset.translations, dataset.n_operations)
    wyckoffs = unsafe_wrap(Vector{Cint}, dataset.wyckoffs, dataset.n_atoms)
    pos = unsafe_wrap(Vector{NTuple{3,Float64}}, dataset.std_positions, dataset.n_std_atoms)
    return Dataset(
        dataset.spacegroup_number,
        dataset.hall_number,
        cchars2string(dataset.international_symbol),
        cchars2string(dataset.hall_symbol),
        cchars2string(dataset.choice),
        tuple2matrix(dataset.transformation_matrix),
        collect(dataset.origin_shift),
        dataset.n_operations,
        rotsFromTuple(r, dataset.n_operations),
        transFromTuple(t, dataset.n_operations),
        dataset.n_atoms,
        [LETTERS[x + 1] for x in wyckoffs],  # Need to add 1 because of C-index starts from 0
        map(
            cchars2string,
            unsafe_wrap(
                Vector{NTuple{7,Cchar}}, dataset.site_symmetry_symbols, dataset.n_atoms
            ),
        ),
        unsafe_wrap(Vector{Cint}, dataset.equivalent_atoms, dataset.n_atoms),
        unsafe_wrap(Vector{Cint}, dataset.crystallographic_orbits, dataset.n_atoms),
        transpose(tuple2matrix(dataset.primitive_lattice)),
        unsafe_wrap(Vector{Cint}, dataset.mapping_to_primitive, dataset.n_atoms),
        dataset.n_std_atoms,
        transpose(tuple2matrix(dataset.std_lattice)),
        unsafe_wrap(Vector{Cint}, dataset.std_types, dataset.n_std_atoms),
        transFromTuple(pos, dataset.n_std_atoms),
        tuple2matrix(dataset.std_rotation_matrix),
        unsafe_wrap(Vector{Cint}, dataset.std_mapping_to_primitive, dataset.n_std_atoms),
        cchars2string(dataset.pointgroup_symbol),
    )
end

function Base.show(io::IO, cell::Cell)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(cell)
        Base.show_default(IOContext(io, :limit => true), cell)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, string(typeof(cell)))
        println(io, " lattice:")
        for row in eachrow(cell.lattice)
            println(io, "  ", join(row, "  "))
        end
        N = natoms(cell)
        println(io, " $N atomic positions:")
        for pos in cell.positions
            println(io, "  ", pos)
        end
        println(io, " $N atoms:")
        println(io, "  ", cell.types)
    end
end
