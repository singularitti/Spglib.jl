using CrystallographyCore: AbstractCell, Cell as CrystallographyCell, basisvectors
using StaticArrays: MVector, SMatrix, SVector
using StructEquality: @struct_hash_equal_isequal

import CrystallographyCore: Lattice, natoms, atomtypes

export Lattice,
    SpglibCell,
    Cell,
    Dataset,
    SpacegroupType,
    basisvectors,
    basis_vectors,
    natoms,
    atomtypes

const basis_vectors = basisvectors  # For backward compatibility

"""
    SpglibCell(lattice, positions, atoms, magmoms=[])

Represent a unit cell with specified lattice, positions, atoms, and magnetic moments.

# Arguments
- `lattice`: lattice of the unit cell. Lattice parameters are given by a ``3×3``
  matrix with floating point values, where ``𝐚``, ``𝐛``, and ``𝐜`` are stored as columns.
  You could also give a vector of 3-vectors, where each vector is a lattice vector.
  See [Basis vectors](@ref) for our conventions and [`Lattice`](@ref) for more examples.
- `positions`: positions of the atoms in the unit cell.
  Fractional atomic positions are given by a vector of ``N`` 3-vectors with floating point
  values, where ``N`` is the number of atoms.
- `atoms`: ``N`` atoms present in the unit cell.
- `magmoms=[]`: magnetic moments on atoms in the unit cell (optional).
  It can be either a vector of ``N`` floating point values for collinear cases or
  a vector of 3-vectors in cartesian coordinates for non-collinear cases.

See also [`Lattice`](@ref).

# Examples
```jldoctest
julia> lattice = Lattice([
           [5.0759761474456697, 5.0759761474456697, 0],
           [-2.8280307701821314, 2.8280307701821314, 0],
           [0, 0, 8.57154746],
       ]);

julia> positions = [
           [0.0, 0.84688439, 0.1203133],
           [0.0, 0.65311561, 0.6203133],
           [0.0, 0.34688439, 0.3796867],
           [0.0, 0.15311561, 0.8796867],
           [0.5, 0.34688439, 0.1203133],
           [0.5, 0.15311561, 0.6203133],
           [0.5, 0.84688439, 0.3796867],
           [0.5, 0.65311561, 0.8796867],
       ];

julia> atoms = fill(35, length(positions));

julia> cell = SpglibCell(lattice, positions, atoms);

julia> lattice = Lattice([
           4 0 0
           0 4 0
           0 0 3
       ]);

julia> positions = [
           [0.0, 0.0, 0.0],
           [0.5, 0.5, 0.5],
           [0.3, 0.3, 0.0],
           [0.7, 0.7, 0.0],
           [0.2, 0.8, 0.5],
           [0.8, 0.2, 0.5],
       ];

julia> atoms = [14, 14, 8, 8, 8, 8];

julia> cell = SpglibCell(lattice, positions, atoms);

julia> lattice = [
           4.0 0.0 0.0
           0.0 4.0 0.0
           0.0 0.0 4.0
       ];

julia> positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];

julia> atoms = [1, 1];

julia> magmoms = [1.0, 1.0];

julia> cell = SpglibCell(lattice, positions, atoms, magmoms);
```
"""
@struct_hash_equal_isequal struct SpglibCell{L,P,T,M} <: AbstractCell
    lattice::Lattice{L}
    positions::Vector{MVector{3,P}}
    atoms::Vector{T}
    magmoms::Vector{M}
end
function SpglibCell(lattice, positions, atoms, magmoms=[])
    if !(lattice isa Lattice)
        lattice = Lattice(lattice)
    end
    N = length(atoms)
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
    L, T, M = eltype(lattice), eltype(atoms), eltype(magmoms)
    return SpglibCell{L,P,T,M}(lattice, positions, atoms, magmoms)
end
SpglibCell(cell::CrystallographyCell, magmoms) =
    SpglibCell(cell.lattice, cell.positions, cell.atoms, magmoms)
const Cell = SpglibCell

"""
    natoms(cell::SpglibCell)

Return the number of atoms in the specified `cell`.
"""
natoms(cell::SpglibCell) = length(cell.atoms)

"""
    atomtypes(cell::SpglibCell)

Return the unique atom types in the specified `cell`.
"""
atomtypes(cell::SpglibCell) = unique(cell.atoms)

"""
    Lattice(cell::SpglibCell)

Get the lattice from a `cell`.
"""
Lattice(cell::SpglibCell) = cell.lattice

# This is an internal function, do not export!
function _expand_cell(cell::SpglibCell)
    lattice, positions, atoms, magmoms = cell.lattice,
    cell.positions, cell.atoms,
    cell.magmoms
    # Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L32-L35 and https://github.com/spglib/spglib/blob/444e061/python/spglib/spglib.py#L953-L975
    lattice = Base.cconvert(Matrix{Cdouble}, transpose(lattice))  # `transpose` must before `cconvert`!
    positions = Base.cconvert(Matrix{Cdouble}, reduce(hcat, positions))
    atomtypes = unique(atoms)
    atoms = collect(Cint, findfirst(==(atom), atomtypes) for atom in atoms)  # Mapping between unique atom types and atom indices
    if !isempty(magmoms)
        magmoms = if eltype(magmoms) <: AbstractVector
            Base.cconvert(Matrix{Cdouble}, reduce(hcat, magmoms))
        else
            Base.cconvert(Vector{Cdouble}, magmoms)
        end
    end
    return lattice, positions, atoms, magmoms
end
function _expand_cell(cell::CrystallographyCell)
    lattice, positions, atoms = cell.lattice, cell.positions, cell.atoms
    lattice = Base.cconvert(Matrix{Cdouble}, transpose(lattice))  # `transpose` must before `cconvert`!
    positions = Base.cconvert(Matrix{Cdouble}, reduce(hcat, positions))
    atomtypes = unique(atoms)
    atoms = collect(Cint, findfirst(==(atom), atomtypes) for atom in atoms)  # Mapping between unique atom types and atom indices
    return lattice, positions, atoms
end

# This is an internal type, do not export!
struct SpglibSpacegroupType
    number::Cint
    international_short::NTuple{11,Cchar}
    international_full::NTuple{20,Cchar}
    international::NTuple{32,Cchar}
    schoenflies::NTuple{7,Cchar}
    hall_number::Cint
    hall_symbol::NTuple{17,Cchar}
    choice::NTuple{6,Cchar}
    pointgroup_international::NTuple{6,Cchar}
    pointgroup_schoenflies::NTuple{4,Cchar}
    arithmetic_crystal_class_number::Cint
    arithmetic_crystal_class_symbol::NTuple{7,Cchar}
end

"""
    SpacegroupType(number, international_short, international_full, international, schoenflies, hall_number, hall_symbol, choice, pointgroup_international, pointgroup_schoenflies, arithmetic_crystal_class_number, arithmetic_crystal_class_symbol)

Represent `SpglibSpacegroupType`, see its [official documentation](https://spglib.github.io/spglib/api.html#spg-get-spacegroup-type).

See also [`get_spacegroup_type`](@ref), [`get_spacegroup_type_from_symmetry`](@ref).
"""
struct SpacegroupType
    number::Int32
    international_short::String
    international_full::String
    international::String
    schoenflies::String
    hall_number::Int32
    hall_symbol::String
    choice::String
    pointgroup_international::String
    pointgroup_schoenflies::String
    arithmetic_crystal_class_number::Int32
    arithmetic_crystal_class_symbol::String
end

# This is an internal type, do not export!
struct SpglibDataset
    spacegroup_number::Cint
    hall_number::Cint
    international_symbol::NTuple{11,Cchar}
    hall_symbol::NTuple{17,Cchar}
    choice::NTuple{6,Cchar}
    transformation_matrix::NTuple{3,NTuple{3,Cdouble}}
    origin_shift::NTuple{3,Cdouble}
    n_operations::Cint
    rotations::Ptr{NTuple{3,NTuple{3,Cint}}}
    translations::Ptr{NTuple{3,Cdouble}}
    n_atoms::Cint
    wyckoffs::Ptr{Cint}
    site_symmetry_symbols::Ptr{NTuple{7,Cchar}}
    equivalent_atoms::Ptr{Cint}
    crystallographic_orbits::Ptr{Cint}
    primitive_lattice::NTuple{3,NTuple{3,Cdouble}}
    mapping_to_primitive::Ptr{Cint}
    n_std_atoms::Cint
    std_lattice::NTuple{3,NTuple{3,Cdouble}}
    std_types::Ptr{Cint}
    std_positions::Ptr{NTuple{3,Cdouble}}
    std_rotation_matrix::NTuple{3,NTuple{3,Cdouble}}
    std_mapping_to_primitive::Ptr{Cint}
    pointgroup_symbol::NTuple{6,Cchar}
end

"""
    Dataset(spacegroup_number, hall_number, international_symbol, hall_symbol, choice, transformation_matrix, origin_shift, n_operations, rotations, translations, n_atoms, wyckoffs, site_symmetry_symbols, equivalent_atoms, crystallographic_orbits, primitive_lattice, mapping_to_primitive, n_std_atoms, std_lattice, std_types, std_positions, std_rotation_matrix, std_mapping_to_primitive, pointgroup_symbol)

Represent `SpglibDataset`, see its [official documentation](https://spglib.github.io/spglib/dataset.html#spglib-dataset).

# Arguments
- `spacegroup_number`: international space group number.
- `hall_number`: Hall number. This number is used in
  [`get_symmetry_from_database`](@ref) and [`get_spacegroup_type`](@ref).
- `international_symbol`: international short symbol.
- `hall_symbol`: Hall symbol.
- `choice`: centring, origin, basis vector setting.
- `transformation_matrix`: See the detail at
  [Transformation matrix and origin shift](@ref).
- `origin shift`: See the detail at
  [Transformation matrix and origin shift](@ref).
- `n_operations`: number of symmetry operations.
- `rotations` and `translations`: rotation matrices and
  translation vectors. See [`get_symmetry`](@ref) for more details.
- `n_atoms`: number of atoms in the input unit cell.
- `wyckoffs`: Wyckoff letters.
- `site_symmetry_symbols`: site-symmetry symbols (experimental).
- `equivalent_atoms`: mapping table to equivalent atoms.
- `crystallographic_orbits` : mapping table to equivalent atoms (see
  [Wyckoff positions and symmetrically equivalent atoms](@ref) for the difference
  between `equivalent_atoms` and `crystallographic_orbits`).
- `primitive_lattice` : basis vectors of a primitive cell.
- `mapping_to_primitive`: mapping table to atoms in the primitive cell.
- `n_std_atoms`: number of atoms in the standardized unit cell.
- `std_lattice`, `std_positions`, `std_types`: standardized
  crystal structure corresponding to the Hall symbol found. These are
  equivalently given in the array formats of `lattice`,
  `positions`, and `atoms` presented at [`SpglibCell`](@ref), respectively.
- `std_rotation_matrix`: see the detail at
  [Standardized crystal structure after idealization](@ref).
- `std_mapping_to_primitive`: Mapping table from atoms in the
  standardized crystal structure to the atoms in the primitive cell.
- `pointgroup_symbol`: symbol of the crystallographic point group in
  the Hermann–Mauguin notation.

See also [`get_dataset`](@ref), [`get_dataset_with_hall_number`](@ref).
"""
@struct_hash_equal_isequal struct Dataset
    spacegroup_number::Int32
    hall_number::Int32
    international_symbol::String
    hall_symbol::String
    choice::String
    transformation_matrix::SMatrix{3,3,Float64,9}
    origin_shift::SVector{3,Float64}
    n_operations::Int32
    rotations::Vector{SMatrix{3,3,Int32,9}}
    translations::Vector{SVector{3,Float64}}
    n_atoms::Int32
    wyckoffs::Vector{Char}
    site_symmetry_symbols::Vector{String}
    equivalent_atoms::Vector{Int32}
    crystallographic_orbits::Vector{Int32}
    primitive_lattice::Lattice{Float64}
    mapping_to_primitive::Vector{Int32}
    n_std_atoms::Int32
    std_lattice::Lattice{Float64}
    std_types::Vector{Int32}
    std_positions::Vector{SVector{3,Float64}}
    std_rotation_matrix::SMatrix{3,3,Float64,9}
    std_mapping_to_primitive::Vector{Int32}
    pointgroup_symbol::String
end

function Base.convert(::Type{Dataset}, dataset::SpglibDataset)
    international_symbol = tostring(dataset.international_symbol)
    hall_symbol = tostring(dataset.hall_symbol)
    choice = tostring(dataset.choice)
    transformation_matrix = transpose(
        _convert(SMatrix{3,3,Float64}, dataset.transformation_matrix)
    )
    rotations = transpose.(
        _convert(SMatrix{3,3,Int32}, unsafe_load(dataset.rotations, i)) for i in Base.OneTo(dataset.n_operations)
    )
    translations = SVector{3}.(
        unsafe_load(dataset.translations, i) for i in Base.OneTo(dataset.n_operations)
    )
    wyckoffs = unsafe_wrap(Vector{Int32}, dataset.wyckoffs, dataset.n_atoms)
    wyckoffs = [('a':'z')[w + 1] for w in wyckoffs]  # Need to add 1 because of C-index starts from 0
    site_symmetry_symbols = tostring.(
        unsafe_load(dataset.site_symmetry_symbols, i) for i in Base.OneTo(dataset.n_atoms)
    )
    equivalent_atoms =  # Need to add 1 because of C-index starts from 0
        unsafe_wrap(Vector{Int32}, dataset.equivalent_atoms, dataset.n_atoms) .+ 1
    crystallographic_orbits =  # Need to add 1 because of C-index starts from 0
        unsafe_wrap(Vector{Int32}, dataset.crystallographic_orbits, dataset.n_atoms) .+ 1
    primitive_lattice = Lattice(
        transpose(_convert(SMatrix{3,3,Float64}, dataset.primitive_lattice))
    )
    mapping_to_primitive =  # Need to add 1 because of C-index starts from 0
        unsafe_wrap(Vector{Int32}, dataset.mapping_to_primitive, dataset.n_atoms) .+ 1
    std_lattice = Lattice(transpose(_convert(SMatrix{3,3,Float64}, dataset.std_lattice)))
    std_types = unsafe_wrap(Vector{Int32}, dataset.std_types, dataset.n_std_atoms)
    std_positions = SVector{3}.(
        unsafe_load(dataset.std_positions, i) for i in Base.OneTo(dataset.n_std_atoms)
    )
    std_rotation_matrix = transpose(
        _convert(SMatrix{3,3,Float64}, dataset.std_rotation_matrix)
    )
    std_mapping_to_primitive =  # Need to add 1 because of C-index starts from 0
        unsafe_wrap(Vector{Int32}, dataset.std_mapping_to_primitive, dataset.n_std_atoms) .+
        1
    pointgroup_symbol = tostring(dataset.pointgroup_symbol)
    return Dataset(
        dataset.spacegroup_number,
        dataset.hall_number,
        international_symbol,
        hall_symbol,
        choice,
        transformation_matrix,
        dataset.origin_shift,
        dataset.n_operations,
        rotations,
        translations,
        dataset.n_atoms,
        wyckoffs,
        site_symmetry_symbols,
        equivalent_atoms,
        crystallographic_orbits,
        primitive_lattice,
        mapping_to_primitive,
        dataset.n_std_atoms,
        std_lattice,
        std_types,
        std_positions,
        std_rotation_matrix,
        std_mapping_to_primitive,
        pointgroup_symbol,
    )
end
function Base.convert(::Type{SpacegroupType}, spgtype::SpglibSpacegroupType)
    international_short = tostring(spgtype.international_short)
    international_full = tostring(spgtype.international_full)
    international = tostring(spgtype.international)
    schoenflies = tostring(spgtype.schoenflies)
    hall_symbol = tostring(spgtype.hall_symbol)
    choice = tostring(spgtype.choice)
    pointgroup_international = tostring(spgtype.pointgroup_international)
    pointgroup_schoenflies = tostring(spgtype.pointgroup_schoenflies)
    arithmetic_crystal_class_symbol = tostring(spgtype.arithmetic_crystal_class_symbol)
    return SpacegroupType(
        spgtype.number,
        international_short,
        international_full,
        international,
        schoenflies,
        spgtype.hall_number,
        hall_symbol,
        choice,
        pointgroup_international,
        pointgroup_schoenflies,
        spgtype.arithmetic_crystal_class_number,
        arithmetic_crystal_class_symbol,
    )
end

_convert(::Type{SMatrix{N,N,T}}, tuple::NTuple{N,NTuple{N,T}}) where {N,T} =
    SMatrix{N,N,T}(Tuple(Iterators.flatten(tuple)))
