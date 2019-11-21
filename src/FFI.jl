"""
# module FFI



# Examples

```jldoctest
julia>
```
"""
module FFI

using Compat: isnothing
using CoordinateTransformations
using Parameters: @unpack
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
       get_multiplicity
    #    get_ir_reciprocal_mesh,
    #    get_stabilized_reciprocal_mesh

# This is an internal type, do not export!
const TupleOrVec = Union{Tuple,AbstractVector}

# This is an internal function, do not export!
function get_ccell(cell::Cell{<:AbstractMatrix,<:AbstractMatrix})
    @unpack lattice, positions, numbers = cell
    # Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L32-L35
    clattice = convert(Matrix{Cdouble}, lattice)
    cpositions = convert(Matrix{Cdouble}, positions)
    cnumbers = Cint[findfirst(isequal(u), unique(numbers)) for u in numbers]
    return Cell(clattice, cpositions, cnumbers)
end
get_ccell(cell::Cell{<:AbstractVector{<:TupleOrVec},<:AbstractVector{<:TupleOrVec}}) = cell

# This is an internal function, do not export!
function trunc_trailing_zeros(s)
    i = findfirst(iszero, s)
    isnothing(i) && return s
    return s[1:findfirst(iszero, s)-1]
end # function trunc_trailing_zeros

# Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L70
# This is an internal function, do not export!
cchars_to_string(s::AbstractVector{Cchar}) =
    convert(Array{Char}, trunc_trailing_zeros(s)) |> join

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
end # function get_symmetry

function get_dataset(cell::Cell; symprec::Real = 1e-8)
    @unpack lattice, positions, numbers = get_ccell(cell)
    # dataset = Wrapper.spg_get_dataset(
    #     lattice,
    #     positions,
    #     numbers,
    #     length(numbers),
    #     symprec,
    # )
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
end # function get_dataset

function get_spacegroup_type(hall_number::Integer)
    spgtype =
        ccall((:spg_get_spacegroup_type, libsymspg), CspaceGroup, (Cint,), hall_number)
    return convert(SpaceGroup, spgtype)
end # function get_spacegroup_type

function get_international(cell::Cell, symprec::Real = 1e-8)
    result = zeros(Cchar, 11)
    @unpack lattice, positions, numbers = get_ccell(cell)
    # exitcode = Wrapper.spg_get_international(
    #     result,
    #     lattice,
    #     positions,
    #     numbers,
    #     length(numbers),
    #     symprec,
    # )
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
    return cchars_to_string(result)
end # function get_international

function get_schoenflies(cell::Cell, symprec::Real = 1e-8)
    result = zeros(Cchar, 11)
    @unpack lattice, positions, numbers = get_ccell(cell)
    # exitcode = Wrapper.spg_get_schoenflies(
    #     result,
    #     lattice,
    #     positions,
    #     numbers,
    #     length(numbers),
    #     symprec,
    # )
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
    return cchars_to_string(result)
end # function get_schoenflies

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
end # function standardize_cell

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
end # function niggli_reduce

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
end # function delaunay_reduce

# function get_ir_reciprocal_mesh(
#     cell::Cell,
#     grid::Vector{T},
#     shift::Vector{T} = [0, 0, 0];
#     is_time_reversal::Bool = true,
#     symprec::Real = 1e-5,
# ) where {T<:Integer}
#     all(
#         x -> x in (zero(T), one(T)),
#         shift,
#     ) || throw(ArgumentError("The shift can be only a vector of ones or zeros!"))

#     qpoints_amount = prod(grid)
#     grid_address = Array{Cint}(undef, qpoints_amount, 3)
#     mapping = Array{Cint}(undef, qpoints_amount)
#     ccell = get_ccell(cell)
#     @unpack lattice, positions, numbers = ccell

#     ret = ccall(
#         (:spg_get_ir_reciprocal_mesh, spglib),
#         Cint,
#         (
#          Ptr{Cint},
#          Ptr{Cint},
#          Ptr{Cint},
#          Ptr{Cint},
#          Cint,
#          Ptr{Cdouble},
#          Ptr{Cdouble},
#          Ptr{Cint},
#          Cint,
#          Cdouble,
#         ),
#         grid_address,
#         mapping,
#         grid,
#         shift,
#         is_time_reversal,
#         lattice,
#         positions,
#         numbers,
#         length(numbers),
#         symprec,
#     )
#     ret != qpoints_amount && error("Something wrong happens when finding mesh!")

#     mapping, grid_address
# end # function get_ir_reciprocal_mesh

# function get_stabilized_reciprocal_mesh(
#     rotations::Vector{Matrix{T}},
#     grid::Vector{T},
#     shift::Vector{T} = [0, 0, 0];
#     qpoints::Vector{} = nothing,
#     is_time_reversal::Bool = true,
# ) where {T<:Integer}
#     all(
#         x -> x in (zero(T), one(T)),
#         shift,
#     ) || throw(ArgumentError("The shift can be only a vector of ones or zeros!"))

#     qpoints_amount = prod(grid)
#     grid_address = Array{Cint}(undef, qpoints_amount, 3)
#     mapping_table = Array{Cint}(undef, qpoints_amount)
#     isnothing(qpoints) ? qpoints = Float64[0, 0, 0] : qpoints = Vector(qpoints)

#     ret = ccall(
#         (:spg_get_stabilized_reciprocal_mesh, spglib),
#         Cint,
#         (
#          Ptr{Cint},
#          Ptr{Cint},
#          Ptr{Cint},
#          Ptr{Cint},
#          Cint,
#          Cint,
#          Ptr{Cint},
#          Cint,
#          Ptr{Cint},
#         ),
#         grid_address,
#         mapping_table,
#         grid,
#         shift,
#         is_time_reversal,
#         length(rotations),
#         rotations,
#         length(qpoints),
#         qpoints,
#     )
#     ret != qpoints_amount && error("Something wrong happens when finding mesh!")

#     mapping_table, grid_address
# end # function get_stabilized_reciprocal_mesh

"""
    get_multiplicity(cell::Cell, symprec = 1e-8)

Return the exact number of symmetry operations. An error is thrown when it failed.
"""
function get_multiplicity(cell::Cell, symprec::Real = 1e-8)
    @unpack lattice, positions, numbers = get_ccell(cell)
    # exitcode = Wrapper.spg_get_international(
    #     result,
    #     lattice,
    #     positions,
    #     numbers,
    #     length(numbers),
    #     symprec,
    # )
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
end # function get_multiplicity

function Base.convert(::Type{T}, dataset::Cdataset) where {T<:Dataset}
    f = name -> getfield(dataset, name) |> convert_field
    return T(map(f, fieldnames(T))...)
end # function Base.convert
function Base.convert(::Type{T}, spgtype::CspaceGroup) where {T<:SpaceGroup}
    f = name -> getfield(spgtype, name) |> convert_field
    # Reference: https://discourse.julialang.org/t/construct-an-immutable-type-from-a-dict/26709/2
    return T(map(f, fieldnames(T))...)
end # function Base.convert

end
