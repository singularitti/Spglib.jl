export get_symmetry,
    get_symmetry_from_database,
    get_dataset,
    get_dataset_with_hall_number,
    get_multiplicity,
    get_international,
    get_schoenflies,
    get_spacegroup_type,
    get_spacegroup_type_from_symmetry,
    get_hall_number_from_symmetry

# See https://github.com/spglib/spglib/blob/444e061/python/spglib/spglib.py#L115-L165
"""
    get_symmetry(cell::AbstractCell, symprec=1e-5)

Return the symmetry operations (rotations, translations) of a `cell`.

Returned value `rotations` is a `Vector` of matrices. It has the length of
"number of symmetry operations". Each matrix is a ``3 \\times 3`` integer matrix.
Returned value `translations` is a `Vector` of vectors. It has the length of
"number of symmetry operations". Each vector is a length-``3`` vector of floating point numbers.

The orders of the rotation matrices and the translation
vectors correspond with each other, e.g., the second symmetry
operation is organized by the set of the second rotation matrix and second
translation vector in the respective arrays. Therefore a set of
symmetry operations may obtained by
`[(r, t) for r, t in zip(rotations, translations)]`.

The operations are given with respect to the fractional coordinates
(not for Cartesian coordinates). The rotation matrix ``\\mathbf{W}`` and translation
vector ``\\text{w}`` are used as follows:

```math
\\tilde{\\mathbf{x}}_{3\\times 1} = \\mathbf{W}_{3\\times 3} \\mathbf{x}_{3\\times 1} + \\text{w}_{3\\times 1}.
```

The three values in the vector are given for the ``a``, ``b``, and ``c`` axes, respectively.

As an exceptional case, if a supercell (or non-primitive cell) has the
basis vectors whose lattice breaks crystallographic point group, the
crystallographic symmetry operations are searched within this broken
symmetry, i.e., at most the crystallographic point group found in this
case is the point group of the lattice. For example, this happens for
the ``2\\times 1\\times 1`` supercell of a conventional cubic unit
cell. This may not be understandable in crystallographic sense, but is
practically useful treatment for research in computational materials
science.

See also [`get_symmetry_with_collinear_spin`](@ref) for magnetic symmetry search.

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

julia> rotations, translations = get_symmetry(cell);

julia> size(rotations) == size(translations) == (16,)
```
"""
function get_symmetry(cell::AbstractCell, symprec=1e-5)
    lattice, positions, atoms = _expand_cell(cell)
    num_atom = natoms(cell)
    # See https://github.com/spglib/spglib/blob/42527b0/python/spglib/spglib.py#L270
    max_size = 48num_atom  # Num of symmetry operations = order of the point group of the space group × num of lattice points
    rotations = Array{Cint,3}(undef, 3, 3, max_size)
    translations = Array{Cdouble,2}(undef, 3, max_size)  # C is row-major order, but Julia is column-major order
    num_sym = @ccall libsymspg.spg_get_symmetry(
        rotations::Ptr{Cint},
        translations::Ptr{Cdouble},
        max_size::Cint,
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        num_atom::Cint,
        symprec::Cdouble,
    )::Cint  # The number of operations is returned.
    check_error()
    rotations, translations = map(
        SMatrix{3,3,Int32,9} ∘ transpose, eachslice(rotations[:, :, 1:num_sym]; dims=3)
    ),  # Remember to transpose, see https://github.com/singularitti/Spglib.jl/blob/8aed6e0/src/core.jl#L195-L198
    map(SVector{3,Float64}, eachcol(translations[:, 1:num_sym]))
    return rotations, translations
end

"""
    get_symmetry_from_database(hall_number)

Return the symmetry operations given a `hall_number`.

This function allows to directly access to the space group operations in the
`spglib` database. To specify the space group type with a specific choice,
`hall_number` is used.

The definition of `hall_number` is found at
[Space group type](https://spglib.readthedocs.io/en/latest/dataset.html#dataset-spg-get-dataset-spacegroup-type).

The serial number from ``1`` to ``530`` which are found at
[list of space groups (Seto's web site)](https://yseto.net/?page_id=29%3E). Be
sure that this is not a standard crystallographic definition.

# Examples
```jldoctest
julia> get_symmetry_from_database(304);
```
"""
function get_symmetry_from_database(hall_number)
    # The maximum number of symmetry operations is 192, see https://github.com/spglib/spglib/blob/77a8e5d/src/spglib.h#L382
    @assert 1 <= hall_number <= 530
    rotations = Array{Cint,3}(undef, 3, 3, 192)
    translations = Array{Cdouble,2}(undef, 3, 192)
    num_sym = @ccall libsymspg.spg_get_symmetry_from_database(
        rotations::Ptr{Cint}, translations::Ptr{Cdouble}, hall_number::Cint
    )::Cint
    check_error()
    rotations, translations = map(
        SMatrix{3,3,Int32,9} ∘ transpose, eachslice(rotations[:, :, 1:num_sym]; dims=3)
    ),  # Remember to transpose, see https://github.com/singularitti/Spglib.jl/blob/8aed6e0/src/core.jl#L195-L198
    map(SVector{3,Float64}, eachcol(translations[:, 1:num_sym]))
    return rotations, translations
end

"""
    get_dataset(cell::AbstractCell, symprec=1e-5)

Search symmetry operations of an input unit cell structure.

For an input unit cell structure, symmetry operations of the crystal
are searched. Then they are compared with the crystallographic
database and the space group type is determined.  The result is
returned as the `Dataset` structure as a dataset.

The detail of the dataset is given at [`Dataset`](@ref).

Dataset corresponding to the space group type in the standard setting
is obtained by `get_dataset`. Here the standard setting means
the first top one among the Hall symbols listed for each space group
type. For example, H setting (hexagonal lattice) is chosen for
rhombohedral crystals. [`get_dataset_with_hall_number`](@ref) explained
below is used to choose different settings such as R setting of
rhombohedral crystals. In this function, the other
crystallographic setting is not obtained.
"""
function get_dataset(cell::AbstractCell, symprec=1e-5)
    lattice, positions, atoms = _expand_cell(cell)
    ptr = @ccall libsymspg.spg_get_dataset(
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        natoms(cell)::Cint,
        symprec::Cdouble,
    )::Ptr{SpglibDataset}
    if ptr == C_NULL
        check_error()
    else
        dataset = unsafe_load(ptr)
        return convert(Dataset, dataset)
    end
end

"""
    get_dataset_with_hall_number(cell::AbstractCell, hall_number::Integer, symprec=1e-5)

Search symmetry operations of an input unit cell structure, using a given Hall number.

For an input unit cell structure, symmetry operations of the crystal
are searched. Then they are compared with the crystallographic
database and the space group type is determined.  The result is
returned as the `Dataset` structure as a dataset.

The detail of the dataset is given at [`Dataset`](@ref).

Dataset corresponding to the space group type in the standard setting
is obtained by `get_dataset`. Here the standard setting means
the first top one among the Hall symbols listed for each space group
type. For example, H setting (hexagonal lattice) is chosen for
rhombohedral crystals. [`get_dataset_with_hall_number`](@ref) explained
below is used to choose different settings such as R setting of
rhombohedral crystals. In this function, the other
crystallographic setting is not obtained.
"""
function get_dataset_with_hall_number(
    cell::AbstractCell, hall_number::Integer, symprec=1e-5
)
    @assert 1 <= hall_number <= 530
    lattice, positions, atoms = _expand_cell(cell)
    ptr = @ccall libsymspg.spg_get_dataset_with_hall_number(
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        natoms(cell)::Cint,
        hall_number::Cint,
        symprec::Cdouble,
    )::Ptr{SpglibDataset}
    if ptr == C_NULL
        check_error()
    else
        dataset = unsafe_load(ptr)
        return convert(Dataset, dataset)
    end
end

"""
    get_multiplicity(cell::AbstractCell, symprec=1e-5)

Return the exact number of symmetry operations.

An error is thrown when it fails.
"""
function get_multiplicity(cell::AbstractCell, symprec=1e-5)
    lattice, positions, atoms = _expand_cell(cell)
    num_sym = @ccall libsymspg.spg_get_multiplicity(
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        natoms(cell)::Cint,
        symprec::Cdouble,
    )::Cint
    check_error()
    return num_sym
end

"""
    get_international(cell::AbstractCell, symprec=1e-5)

Return the space group type in Hermann–Mauguin (international) notation.
"""
function get_international(cell::AbstractCell, symprec=1e-5)
    lattice, positions, atoms = _expand_cell(cell)
    symbol = Vector{Cchar}(undef, 11)
    @ccall libsymspg.spg_get_international(
        symbol::Ptr{Cchar},
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        natoms(cell)::Cint,
        symprec::Cdouble,
    )::Cint
    check_error()
    return tostring(symbol)
end

"""
    get_schoenflies(cell::AbstractCell, symprec=1e-5)

Return the space group type in Schoenflies notation.
"""
function get_schoenflies(cell::AbstractCell, symprec=1e-5)
    lattice, positions, atoms = _expand_cell(cell)
    symbol = Vector{Cchar}(undef, 7)
    @ccall libsymspg.spg_get_schoenflies(
        symbol::Ptr{Cchar},
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        natoms(cell)::Cint,
        symprec::Cdouble,
    )::Cint
    check_error()
    return tostring(symbol)
end

"""
    get_spacegroup_type(hall_number)

Translate Hall number to space group type information.
"""
function get_spacegroup_type(hall_number)
    @assert 1 <= hall_number <= 530
    spgtype = @ccall libsymspg.spg_get_spacegroup_type(
        hall_number::Cint
    )::SpglibSpacegroupType
    return convert(SpacegroupType, spgtype)
end

"""
    get_spacegroup_type_from_symmetry(cell::AbstractCell, symprec=1e-5)

Return space-group type information from symmetry operations.

This is the replacement of [`get_hall_number_from_symmetry`](@ref).

This is expected to work well for the set of symmetry operations whose
distortion is small. The aim of making this feature is to find
space-group-type for the set of symmetry operations given by the other
source than Spglib.

The `SpacegroupType` structure is explained at [`SpacegroupType`](@ref).
The parameter `lattice` is used as the distance measure for `symprec`. If it
is unknown, the following may be a reasonable choice:

```jldoctest
julia> lattice = Lattice([
    1 0 0
    0 1 0
    0 0 1
]);
```
"""
function get_spacegroup_type_from_symmetry(cell::AbstractCell, symprec=1e-5)
    rotations, translations = get_symmetry(cell, symprec)
    num_sym = length(translations)
    rotations, translations = cat(transpose.(rotations)...; dims=3),
    reduce(hcat, translations)
    lattice, _, _ = _expand_cell(cell)
    spgtype = @ccall libsymspg.spg_get_spacegroup_type_from_symmetry(
        rotations::Ptr{Cint},
        translations::Ptr{Cdouble},
        num_sym::Cint,
        lattice::Ptr{Cdouble},
        symprec::Cdouble,
    )::SpglibSpacegroupType
    return convert(SpacegroupType, spgtype)
end

"""
    get_hall_number_from_symmetry(cell::AbstractCell, symprec=1e-5)

Obtain `hall_number` from the set of symmetry operations.

This is expected to work well for the set of symmetry operations whose
distortion is small. The aim of making this feature is to find
space-group-type for the set of symmetry operations given by the other
source than Spglib.

Note that the definition of `symprec` is different from usual one, but is given in the
fractional coordinates and so it should be small like `1e-5`.

!!! warning
    This function will be replaced by [`get_spacegroup_type_from_symmetry`](@ref).
"""
function get_hall_number_from_symmetry(cell::AbstractCell, symprec=1e-5)
    rotations, translations = get_symmetry(cell, symprec)
    num_sym = length(translations)
    rotations, translations = cat(transpose.(rotations)...; dims=3),
    reduce(hcat, translations)
    hall_number = @ccall libsymspg.spg_get_hall_number_from_symmetry(
        rotations::Ptr{Cint}, translations::Ptr{Cdouble}, num_sym::Cint, symprec::Cdouble
    )::Cint
    check_error()
    return hall_number
end

@deprecate get_hall_number_from_symmetry get_spacegroup_type_from_symmetry
