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

import ..Cell
using ..Wrapper

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
       get_ir_reciprocal_mesh,
       get_stabilized_reciprocal_mesh

const TupleOrVec = Union{Tuple,AbstractVector}

function get_ccell(cell::Cell{<:AbstractMatrix,<:AbstractMatrix})
    @unpack lattice, positions, numbers = cell
    # Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L32-L35
    clattice = Iterators.partition(Cdouble.(lattice'), 3) |> collect
    cpositions = Iterators.partition(Cdouble.(positions'), 3) |> collect
    cnumbers = Cint[findfirst(isequal(u), unique(numbers)) for u in numbers]
    return Cell(clattice, cpositions, cnumbers)
end
get_ccell(cell::Cell{<:AbstractVector{<:TupleOrVec},<:AbstractVector{<:TupleOrVec}}) = cell

# Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L70
cchars_to_string(s::AbstractVector{Cchar}) =
    convert(Array{Char}, s[1:findfirst(iszero, s)-1]) |> join

function get_symmetry(cell::Cell, symprec::Real = 1e-8)
    # ccell = get_ccell(cell)
    @unpack lattice, positions, numbers = cell

    maxsize = 48 * length(positions)
    rotations = Array{Int32}(undef, maxsize, 3, 3)
    translations = Array{Float64}(undef, maxsize, 3)

    # numops = ccall(
    #     (:spg_get_symmetry, spglib),
    #     Cint,
    #     (
    #      Ptr{Cint},
    #      Ptr{Cdouble},
    #      Cint,
    #      Ptr{Cdouble},
    #      Ptr{Cdouble},
    #      Ptr{Cint},
    #      Cint,
    #      Cdouble
    #     ),
    #     rotations,
    #     translations,
    #     maxsize,
    #     lattice,
    #     positions,
    #     numbers,
    #     length(numbers),
    #     symprec
    # )
    numops = Wrapper.spg_get_symmetry(
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

# function get_dataset(cell::Cell; symprec::Real = 1e-8)
#     ccell = get_ccell(cell)
#     @unpack lattice, positions, numbers = ccell

#     dataset = ccall(
#         (:spg_get_dataset, spglib),
#         Dataset,
#         (
#          Ptr{Cdouble},
#          Ptr{Cdouble},
#          Ptr{Cint},
#          Cint,
#          Cdouble
#         ),
#         lattice,
#         positions,
#         numbers,
#         length(numbers),
#         symprec
#     )

#     return dataset
# end # function get_dataset

# function get_spacegroup_type(hall_number::Integer)
#     spgtype = ccall(
#         (:spg_get_spacegroup_type, spglib), SpacegroupType, (Cint,), Int32(hall_number)
#     )
#     return spgtype
# end # function get_spacegroup_type

function get_international(cell::Cell, symprec::Real = 1e-8)
    result = zeros(Cchar, 11)
    @unpack lattice, positions, numbers = get_ccell(cell)
    numops = Wrapper.spg_get_international(
        result,
        lattice,
        positions,
        numbers,
        length(numbers),
        symprec,
    )
    numops == 0 && error("Could not determine the international symbol!")
    return cchars_to_string(result)
end # function get_international

function get_schoenflies(cell::Cell, symprec::Real = 1e-8)
    result = zeros(Cchar, 11)
    @unpack lattice, positions, numbers = get_ccell(cell)
    numops = Wrapper.spg_get_schoenflies(
        result,
        lattice,
        positions,
        numbers,
        length(numbers),
        symprec,
    )
    numops == 0 && error("Could not determine the Schoenflies symbol!")
    return cchars_to_string(result)
end # function get_schoenflies

function standardize_cell(
    cell::Cell;
    to_primitive::Bool = false,
    no_idealize::Bool = false,
    symprec::Real = 1e-5,
)
    ccell = get_ccell(cell)
    @unpack lattice, positions, numbers = ccell
    to_primitive, no_idealize = map(x -> convert(Cint, x), (to_primitive, no_idealize))

    atoms_amount = ccall(
        (:spg_standardize_cell, spglib),
        Cint,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Cint, Cdouble),
        lattice,
        positions,
        numbers,
        length(numbers),
        to_primitive,
        no_idealize,
        symprec,
    )
    atoms_amount == 0 && error("Standardizing cell failed!")

    Cell(lattice, positions, numbers)
end # function standardize_cell

find_primitive(cell::Cell, symprec::Real = 1e-5) =
    standardize_cell(cell; to_primitive = true, no_idealize = false, symprec = symprec)

refine_cell(cell::Cell, symprec::Real = 1e-5) =
    standardize_cell(cell; to_primitive = false, no_idealize = false, symprec = symprec)

function niggli_reduce(cell::Cell, symprec::Real = 1e-5)
    # Equivalent to `np.transpose` in https://github.com/atztogo/spglib/blob/f8ddf5b/python/spglib/spglib.py#L869
    lattice = Iterators.partition(getfield(cell, :lattice), 3) |> collect
    # Make sure the `symprec` is a float
    code = Wrapper.niggli_reduce(lattice, float(symprec))  # The result is reassigned to `lattice`.
    iszero(code) && error("Niggli reduce failed!")
    # Equivalent to `np.transpose` in https://github.com/atztogo/spglib/blob/f8ddf5b/python/spglib/spglib.py#L876
    return @set cell.lattice =
        Iterators.flatten(lattice) |> collect |> x -> reshape(x, 3, 3)
end # function niggli_reduce

function delaunay_reduce(cell::Cell, symprec::Real = 1e-5)
    # Equivalent to `np.transpose` in https://github.com/atztogo/spglib/blob/f8ddf5b/python/spglib/spglib.py#L832
    lattice = Iterators.partition(getfield(cell, :lattice), 3) |> collect
    # Make sure the `symprec` is a float
    code = Wrapper.niggli_reduce(lattice, float(symprec))  # The result is reassigned to `lattice`.
    iszero(code) && error("Delaunay reduce failed!")
    # Equivalent to `np.transpose` in https://github.com/atztogo/spglib/blob/f8ddf5b/python/spglib/spglib.py#L840
    return @set cell.lattice =
        Iterators.flatten(lattice) |> collect |> x -> reshape(x, 3, 3)
end # function delaunay_reduce

function get_ir_reciprocal_mesh(
    cell::Cell,
    grid::Vector{T},
    shift::Vector{T} = [0, 0, 0];
    is_time_reversal::Bool = true,
    symprec::Real = 1e-5,
) where {T<:Integer}
    all(
        x -> x in (zero(T), one(T)),
        shift,
    ) || throw(ArgumentError("The shift can be only a vector of ones or zeros!"))

    qpoints_amount = prod(grid)
    grid_address = Array{Cint}(undef, qpoints_amount, 3)
    mapping = Array{Cint}(undef, qpoints_amount)
    ccell = get_ccell(cell)
    @unpack lattice, positions, numbers = ccell

    ret = ccall(
        (:spg_get_ir_reciprocal_mesh, spglib),
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
    ret != qpoints_amount && error("Something wrong happens when finding mesh!")

    mapping, grid_address
end # function get_ir_reciprocal_mesh

function get_stabilized_reciprocal_mesh(
    rotations::Vector{Matrix{T}},
    grid::Vector{T},
    shift::Vector{T} = [0, 0, 0];
    qpoints::Vector{} = nothing,
    is_time_reversal::Bool = true,
) where {T<:Integer}
    all(
        x -> x in (zero(T), one(T)),
        shift,
    ) || throw(ArgumentError("The shift can be only a vector of ones or zeros!"))

    qpoints_amount = prod(grid)
    grid_address = Array{Cint}(undef, qpoints_amount, 3)
    mapping_table = Array{Cint}(undef, qpoints_amount)
    isnothing(qpoints) ? qpoints = Float64[0, 0, 0] : qpoints = Vector(qpoints)

    ret = ccall(
        (:spg_get_stabilized_reciprocal_mesh, spglib),
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
        mapping_table,
        grid,
        shift,
        is_time_reversal,
        length(rotations),
        rotations,
        length(qpoints),
        qpoints,
    )
    ret != qpoints_amount && error("Something wrong happens when finding mesh!")

    mapping_table, grid_address
end # function get_stabilized_reciprocal_mesh

end
