# Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L70
# This is an internal function, do not export!
function cchars2string(itr)
    vec = collect(Char, Iterators.filter(!iszero, itr))
    return String(vec)
end

# See https://github.com/spglib/spglib/blob/444e061/python/spglib/spglib.py#L115-L165
function get_symmetry(cell::Cell, symprec = 1e-5)
    max_size = length(cell.types) * 48
    rotation = Array{Cint,3}(undef, 3, 3, max_size)
    translation = Array{Cdouble,2}(undef, 3, max_size)
    if iszero(cell.magmoms)
        numops = get_symmetry!(rotation, translation, max_size, cell, symprec)
    else
        equivalent_atoms = zeros(length(cell.magmoms))
        primitive_lattice = zeros(Cdouble, 3, 3)
        if ndims(cell.magmoms) == 1
            spin_flips = zeros(length(rotation))
        else
            spin_flips = nothing
        end
        # TODO: unfinished!
    end
    return rotation[:, :, 1:numops], translation[:, 1:numops]
end

function get_symmetry!(
    rotation::AbstractArray,
    translation::AbstractMatrix,
    cell::Cell,
    symprec = 1e-5,
)
    if size(rotation, 3) != size(translation, 2)
        throw(DimensionMismatch("`rotation` & `translation` have different max size!"))
    end
    if !(size(rotation, 1) == size(rotation, 2) == size(translation, 1) == 3)
        throw(ArgumentError("`rotation` & `translation` don't have the right size"))
    end
    @unpack lattice, positions, types = _expand_cell(cell)
    rotation = Base.cconvert(Array{Cint,3}, rotation)
    translation = Base.cconvert(Matrix{Cdouble}, translation)
    max_size = Base.cconvert(Cint, size(rotation, 3))
    number = Base.cconvert(Cint, length(types))
    num_sym = ccall(
        (:spg_get_symmetry, libsymspg),
        Cint,
        (
            Ptr{Cint},
            Ptr{Float64},
            Cint,
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Cint},
            Cint,
            Float64,
        ),
        rotation,
        translation,
        max_size,
        lattice,
        positions,
        types,
        number,
        symprec,
    )
    num_sym == 0 && error("`spg_get_symmetry` failed!")
    return rotation[:, :, 1:num_sym], translation[:, 1:num_sym]
end

function get_symmetry_with_collinear_spin!(
    rotation::AbstractArray{T,3},
    translation::AbstractMatrix,
    equivalent_atoms::AbstractVector,
    max_size::Integer,
    cell::Cell,
    symprec = 1e-5,
) where {T}
    @unpack lattice, positions, types, magmoms = _expand_cell(cell)
    rotation = Base.cconvert(Array{Cint,3}, rotation)
    translation = Base.cconvert(Matrix{Cdouble}, translation)
    equivalent_atoms = Base.cconvert(Vector{Cint}, equivalent_atoms)
    max_size = Base.cconvert(Cint, max_size)
    number = Base.cconvert(Cint, length(types))
    num_sym = ccall(
        (:spg_get_symmetry_with_collinear_spin, libsymspg),
        Cint,
        (
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cint},
            Cint,
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cint},
            Ptr{Cdouble},
            Cint,
            Cdouble,
        ),
        rotation,
        translation,
        equivalent_atoms,
        max_size,
        lattice,
        positions,
        types,
        magmoms,
        number,
        symprec,
    )
    num_sym == 0 && error("`spg_get_symmetry` failed!")
    return num_sym
end
function get_symmetry_with_collinear_spin(cell::Cell, symprec = 1e-5)
    number = length(cell.types)
    max_size = number * 48
    rotation = Array{Cint,3}(undef, 3, 3, max_size)
    translation = Matrix{Cdouble}(undef, 3, max_size)
    equivalent_atoms = Vector{Cint}(undef, number)
    num_sym = get_symmetry_with_collinear_spin!(
        rotation,
        translation,
        equivalent_atoms,
        max_size,
        cell,
        symprec,
    )
    return rotation[:, :, 1:num_sym], translation[:, 1:num_sym], equivalent_atoms
end

"""
    get_hall_number_from_symmetry(rotation::AbstractArray{T,3}, translation::AbstractMatrix, num_operations::Integer, symprec=1e-5) where {T}

Obtain `hall_number` from the set of symmetry operations.

This is expected to work well for the set of symmetry operations whose
distortion is small. The aim of making this feature is to find space-group-type
for the set of symmetry operations given by the other source than spglib. Note
that the definition of `symprec` is different from usual one, but is given in the
fractional coordinates and so it should be small like `1e-5`.
"""
function get_hall_number_from_symmetry(
    rotation::AbstractArray{T,3},
    translation::AbstractMatrix,
    num_operations::Integer,
    symprec = 1e-5,
) where {T}
    rotation = Base.cconvert(Array{Cint,3}, rotation)
    translation = Base.cconvert(Matrix{Cdouble}, translation)
    num_operations = Base.cconvert(Cint, num_operations)
    return ccall(
        (:spg_get_hall_number_from_symmetry, libsymspg),
        Cint,
        (Ptr{Cint}, Ptr{Float64}, Cint, Float64),
        rotation,
        translation,
        num_operations,
        symprec,
    )
end

"""
    get_multiplicity(cell::Cell, symprec=1e-5)

Return the exact number of symmetry operations. An error is thrown when it fails.
"""
function get_multiplicity(cell::Cell, symprec = 1e-5)
    @unpack lattice, positions, types = _expand_cell(cell)
    number = Base.cconvert(Cint, length(types))
    nsymops = ccall(
        (:spg_get_multiplicity, libsymspg),
        Cint,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        lattice,
        positions,
        types,
        number,
        symprec,
    )
    nsymops == 0 && error("Could not determine the multiplicity!")
    return nsymops
end

"""
    get_dataset(cell::Cell, symprec=1e-5)

Search symmetry operations of an input unit cell structure.
"""
function get_dataset(cell::Cell, symprec = 1e-5)
    @unpack lattice, positions, types = _expand_cell(cell)
    number = Base.cconvert(Cint, length(types))
    ptr = ccall(
        (:spg_get_dataset, libsymspg),
        Ptr{SpglibDataset},
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        lattice,
        positions,
        types,
        number,
        symprec,
    )
    raw = unsafe_load(ptr)
    return convert(Dataset, raw)
end

"""
    get_spacegroup_type(hall_number::Integer)

Translate Hall number to space group type information.
"""
function get_spacegroup_type(hall_number::Integer)
    spgtype = ccall(
        (:spg_get_spacegroup_type, libsymspg),
        SpglibSpacegroupType,
        (Cint,),
        hall_number,
    )
    return convert(SpacegroupType, spgtype)
end
"""
    get_spacegroup_type(cell::Cell, symprec=1e-5)

Get `SpacegroupType` from `cell`.
"""
function get_spacegroup_type(cell::Cell, symprec = 1e-5)  # See https://github.com/spglib/spglib/blob/444e061/python/spglib/spglib.py#L307-L324
    dataset = get_dataset(cell, symprec)
    return get_spacegroup_type(dataset.hall_number)
end

"""
    get_international(cell::Cell, symprec=1e-5)

Return the space group type in Hermann–Mauguin (international) notation.
"""
function get_international(cell::Cell, symprec = 1e-5)
    @unpack lattice, positions, types = _expand_cell(cell)
    symbol = Vector{Cchar}(undef, 11)
    exitcode = ccall(
        (:spg_get_international, libsymspg),
        Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        symbol,
        lattice,
        positions,
        types,
        length(types),
        symprec,
    )
    exitcode == 0 && error("Could not determine the international symbol!")
    return cchars2string(symbol)
end

"""
    get_schoenflies(cell::Cell, symprec=1e-5)

Return the space group type in Schoenflies notation.
"""
function get_schoenflies(cell::Cell, symprec = 1e-5)
    @unpack lattice, positions, types = _expand_cell(cell)
    symbol = Vector{Cchar}(undef, 7)
    exitcode = ccall(
        (:spg_get_schoenflies, libsymspg),
        Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        symbol,
        lattice,
        positions,
        types,
        length(types),
        symprec,
    )
    exitcode == 0 && error("Could not determine the Schoenflies symbol!")
    return cchars2string(symbol)
end

function Base.convert(::Type{SpacegroupType}, spgtype::SpglibSpacegroupType)
    values = map(fieldnames(SpacegroupType)) do name
        value = getfield(spgtype, name)
        if value isa Cint
            value
        elseif value isa NTuple{N,Cchar} where {N}
            cchars2string(value)
        else  # This should never happen!
            error("unexpected field type $(typeof(value))!")
        end
    end
    return SpacegroupType(values...)
end
