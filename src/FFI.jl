module FFI

using Compat: isnothing
using CoordinateTransformations
using UnPack: @unpack
using Setfield: @set
using spglib_jll: libsymspg

using ..DataModel: Cell, Dataset, Cdataset, SpaceGroup, CspaceGroup

export get_symmetry,
    get_dataset,
    get_spacegroup_type,
    get_international,
    get_schoenflies,
    standardize_cell,
    find_primitive,
    refine_cell,
    niggli_reduce,
    delaunay_reduce,
    get_multiplicity,
    get_ir_reciprocal_mesh,
    get_stabilized_reciprocal_mesh

# This is an internal function, do not export!
function get_ccell(cell::Cell{<:AbstractMatrix,<:AbstractMatrix})
    @unpack lattice, positions, numbers = cell
    # Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L32-L35
    clattice = convert(Matrix{Cdouble}, lattice)
    cpositions = convert(Matrix{Cdouble}, positions)
    cnumbers = Cint[findfirst(isequal(u), unique(numbers)) for u in numbers]
    return Cell(clattice, cpositions, cnumbers)
end

# This is an internal function, do not export!
trunc_trailing_zeros(vec) = Iterators.filter(!iszero, vec)

# Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L70
# This is an internal function, do not export!
cchars2string(vec::AbstractVector{Cchar}) =
    String(convert(Vector{Char}, trunc_trailing_zeros(vec)))

convert_field(x) = x  # Integers
convert_field(x::NTuple{N,Integer}) where {N} = String(collect(trunc_trailing_zeros(x)))
convert_field(x::Ptr) = unsafe_load(x)
convert_field(x::NTuple{M,NTuple{N,Number}}) where {M,N} =
    transpose(reshape(collect(Iterators.flatten(x)), N, M))
convert_field(x::NTuple{N,Number}) where {N} = collect(x)

function get_symmetry(cell::Cell, symprec::Real = 1e-8)
    @unpack lattice, positions, numbers = get_ccell(cell)
    maxsize = 52 * length(positions)
    rotations = Array{Cint}(undef, maxsize, 3, 3)
    translations = Array{Cdouble}(undef, maxsize, 3)
    numops = ccall(
        (:spg_get_symmetry, libsymspg),
        Cint,
        (
            Ptr{Cint},
            Ptr{Cdouble},
            Cint,
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cint},
            Cint,
            Cdouble,
        ),
        rotations,
        translations,
        maxsize,
        lattice,
        positions,
        numbers,
        length(numbers),
        symprec,
    )
    numops == 0 && error("Could not determine symmetries!")
    [AffineMap(transpose(rotations[:, :, i]), translations[:, i]) for i in 1:numops]
end

function get_dataset(cell::Cell; symprec::Real = 1e-8)
    @unpack lattice, positions, numbers = get_ccell(cell)
    ptr = ccall(
        (:spg_get_dataset, libsymspg),
        Ptr{Cdataset},
        (Ptr{NTuple{3,Cdouble}}, Ptr{NTuple{3,Cdouble}}, Ptr{Cint}, Cint, Cdouble),
        lattice,
        positions,
        numbers,
        length(numbers),
        symprec,
    )
    return convert(Dataset, unsafe_load(ptr))
end

function get_spacegroup_type(hall_number::Integer)
    spgtype =
        ccall((:spg_get_spacegroup_type, libsymspg), CspaceGroup, (Cint,), hall_number)
    return convert(SpaceGroup, spgtype)
end

function get_international(cell::Cell, symprec::Real = 1e-8)
    result = zeros(Cchar, 11)
    @unpack lattice, positions, numbers = get_ccell(cell)
    exitcode = ccall(
        (:spg_get_international, libsymspg),
        Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        result,
        lattice,
        positions,
        numbers,
        length(numbers),
        symprec,
    )
    exitcode == 0 && error("Could not determine the international symbol!")
    return cchars2string(result)
end

function get_schoenflies(cell::Cell, symprec::Real = 1e-8)
    result = zeros(Cchar, 11)
    @unpack lattice, positions, numbers = get_ccell(cell)
    exitcode = ccall(
        (:spg_get_schoenflies, libsymspg),
        Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        result,
        lattice,
        positions,
        numbers,
        length(numbers),
        symprec,
    )
    exitcode == 0 && error("Could not determine the Schoenflies symbol!")
    return cchars2string(result)
end

function standardize_cell(
    cell::Cell;
    to_primitive::Bool = false,
    no_idealize::Bool = false,
    symprec::Real = 1e-5,
)
    @unpack lattice, positions, numbers = get_ccell(cell)
    exitcode = ccall(
        (:spg_standardize_cell, libsymspg),
        Cint,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Cint, Cdouble),
        lattice,
        position,
        numbers,
        length(numbers),
        to_primitive,
        no_idealize,
        symprec,
    )
    exitcode == 0 && error("Standardizing cell failed!")
    return Cell(lattice, positions, numbers)  # They have been changed now.
end

find_primitive(cell::Cell, symprec::Real = 1e-5) =
    standardize_cell(cell; to_primitive = true, no_idealize = false, symprec = symprec)

refine_cell(cell::Cell, symprec::Real = 1e-5) =
    standardize_cell(cell; to_primitive = false, no_idealize = false, symprec = symprec)

function niggli_reduce(cell::Cell, symprec::Real = 1e-5)
    # Equivalent to `np.transpose` in https://github.com/atztogo/spglib/blob/f8ddf5b/python/spglib/spglib.py#L869
    lattice = Iterators.partition(getfield(cell, :lattice), 3) |> collect
    # The result is reassigned to `lattice`.
    exitcode = ccall(
        (:spg_niggli_reduce, libsymspg),
        Cint,
        (Ptr{Cdouble}, Cdouble),
        lattice,
        symprec,
    )
    iszero(exitcode) && error("Niggli reduce failed!")
    # Equivalent to `np.transpose` in https://github.com/atztogo/spglib/blob/f8ddf5b/python/spglib/spglib.py#L876
    return @set cell.lattice =
        Iterators.flatten(lattice) |> collect |> x -> reshape(x, 3, 3)
end

function delaunay_reduce(cell::Cell, symprec::Real = 1e-5)
    # Equivalent to `np.transpose` in https://github.com/atztogo/spglib/blob/f8ddf5b/python/spglib/spglib.py#L832
    lattice = Iterators.partition(getfield(cell, :lattice), 3) |> collect
    # The result is reassigned to `lattice`.
    exitcode = ccall(
        (:spg_delaunay_reduce, libsymspg),
        Cint,
        (Ptr{Cdouble}, Cdouble),
        lattice,
        symprec,
    )
    iszero(exitcode) && error("Delaunay reduce failed!")
    # Equivalent to `np.transpose` in https://github.com/atztogo/spglib/blob/f8ddf5b/python/spglib/spglib.py#L840
    return @set cell.lattice =
        Iterators.flatten(lattice) |> collect |> x -> reshape(x, 3, 3)
end

function get_ir_reciprocal_mesh(
    cell::Cell,
    grid::AbstractVector{<:Integer},
    shift::AbstractVector{<:Integer} = [0, 0, 0];
    is_time_reversal::Bool = true,
    symprec::Real = 1e-5,
)
    @assert(length(grid) == length(shift) == 3)
    @assert(all(isone(x) || iszero(x) for x in shift))
    npoints = prod(grid)
    grid_address = Matrix{Cint}(undef, npoints, 3)
    mapping = Vector{Cint}(undef, npoints)
    @unpack lattice, positions, numbers = get_ccell(cell)
    exitcode = ccall(
        (:spg_get_ir_reciprocal_mesh, libsymspg),
        Cint,
        (
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Cint,
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cint},
            Cint,
            Cdouble,
        ),
        grid_address,
        mapping,
        grid,
        shift,
        is_time_reversal,
        lattice,
        positions,
        numbers,
        length(numbers),
        symprec,
    )
    @assert(exitcode > 0, "Something wrong happens when finding mesh!")
    return mapping, grid_address
end

function get_stabilized_reciprocal_mesh(
    rotations::AbstractVector{AbstractMatrix{<:Integer}},
    grid::AbstractVector{<:Integer},
    shift::AbstractVector{<:Integer} = [0, 0, 0];
    qpoints::AbstractMatrix{<:AbstractFloat} = [[0, 0, 0]],
    is_time_reversal::Bool = true,
)
    @assert(length(grid) == length(shift) == 3)
    @assert(all(isone(x) || iszero(x) for x in shift))
    @assert(size(qpoints, 2) == 3)
    npoints = prod(grid)
    grid_address = Matrix{Cint}(undef, npoints, 3)
    mapping = Vector{Cint}(undef, npoints)
    exitcode = ccall(
        (:spg_get_stabilized_reciprocal_mesh, libsymspg),
        Cint,
        (
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Cint,
            Cint,
            Ptr{Cint},
            Cint,
            Ptr{Cint},
        ),
        grid_address,
        mapping,
        grid,
        shift,
        is_time_reversal,
        length(rotations),
        rotations,
        length(qpoints),
        qpoints,
    )
    @assert(exitcode > 0, "Something wrong happens when finding mesh!")
    return mapping, grid_address
end

"""
    get_multiplicity(cell::Cell, symprec = 1e-8)

Return the exact number of symmetry operations. An error is thrown when it fails.
"""
function get_multiplicity(cell::Cell, symprec::Real = 1e-8)
    @unpack lattice, positions, numbers = get_ccell(cell)
    nsymops = ccall(
        (:spg_get_multiplicity, libsymspg),
        Cint,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        lattice,
        positions,
        numbers,
        length(numbers),
        symprec,
    )
    nsymops == 0 && error("Could not determine the multiplicity!")
    return nsymops
end

function Base.convert(::Type{T}, dataset::Cdataset) where {T<:Dataset}
    f = name -> getfield(dataset, name) |> convert_field
    return T(map(f, fieldnames(T))...)
end
function Base.convert(::Type{T}, spgtype::CspaceGroup) where {T<:SpaceGroup}
    f = name -> getfield(spgtype, name) |> convert_field
    # Reference: https://discourse.julialang.org/t/construct-an-immutable-type-from-a-dict/26709/2
    return T(map(f, fieldnames(T))...)
end

end
