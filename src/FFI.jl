"""
# module FFI



# Examples

```jldoctest
julia>
```
"""
module FFI

using CoordinateTransformations
using DataStructures: counter
using Setfield: @set

using SpgLib.DataModel: Cell

export get_symmetry, get_international, get_schoenflies,
    standardize_cell, find_primitive, refine_cell,
    niggli_reduce, delaunay_reduce

include(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))

const LIBVERSION = VersionNumber(
  ccall((:spg_get_major_version, spglib), Cint, ()),
  ccall((:spg_get_minor_version, spglib), Cint, ()),
  ccall((:spg_get_micro_version, spglib), Cint, ()),
)

function getfields(obj, fields...)
    Tuple(getfield(obj, name) for name in fields)
end

function get_ccell(cell::Cell)::Cell
    lattice, positions, types = getfields(cell, :lattice, :positions, :numbers)
    clattice = convert(Matrix{Cdouble}, lattice)
    cpositions = convert(Matrix{Cdouble}, positions)
    cnumbers = convert(Vector{Cint}, [repeat([i], v) for (i, v) in (enumerate ∘ values ∘ counter)(types)] |> Iterators.flatten |> collect)
    return Cell(clattice, cpositions, cnumbers)
end

cchars_to_string(s::Vector{Cchar}) = map(Char, s) |> join |> x -> split(x, "\0") |> first

function get_symmetry(cell::Cell; symprec::Real = 1e-8)
    size(cell.positions, 2) != length(cell.numbers) && throw(DimensionMismatch("The number of positions and atomic types do not match!"))
    size(cell.positions, 1) != 3 && error("Operations in 3D space is supported here!")

    maxsize = 52
    rotations = Array{Cint}(undef, 3, 3, maxsize)
    translations = Array{Cdouble}(undef, 3, maxsize)

    ccell = get_ccell(cell)
    clattice, cpositions, cnumbers = getfields(ccell, :lattice, :positions, :numbers)

    numops = ccall((:spg_get_symmetry, spglib), Cint,
        (Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        rotations, translations, maxsize, clattice, cpositions, cnumbers, length(cnumbers), symprec)
    numops == 0 && error("Could not determine symmetries!")

    [AffineMap(transpose(rotations[:, :, i]), translations[:, i]) for i in 1:numops]
end

function get_international(cell::Cell; symprec::Real = 1e-8)
    result = zeros(Cchar, 11)

    ccell = get_ccell(cell)
    clattice, cpositions, cnumbers = getfields(ccell, :lattice, :positions, :numbers)

    numops = ccall((:spg_get_international, spglib), Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        result, clattice, cpositions, cnumbers, length(cnumbers), symprec)
    numops == 0 && error("Could not determine the international symbol!")

    cchars_to_string(result)
end

function get_schoenflies(cell::Cell; symprec::Real = 1e-8)
    result = zeros(Cchar, 11)

    ccell = get_ccell(cell)
    clattice, cpositions, cnumbers = getfields(ccell, :lattice, :positions, :numbers)

    numops = ccall((:spg_get_schoenflies, spglib), Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        result, clattice, cpositions, cnumbers, length(cnumbers), symprec)
    numops == 0 && error("Could not determine the Schoenflies symbol!")

    cchars_to_string(result)
end

function standardize_cell(cell::Cell, to_primitive::Bool = false, no_idealize::Bool = false, symprec::Real = 1e-5)
    ccell = get_ccell(cell)
    clattice, cpositions, cnumbers = getfields(ccell, :lattice, :positions, :numbers)
    to_primitive, no_idealize = map(x -> convert(Cint, x), (to_primitive, no_idealize))

    atoms_amount = ccall((:spg_standardize_cell, spglib), Cint,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Cint, Cdouble),
        clattice, cpositions, cnumbers, length(cnumbers), to_primitive, no_idealize, symprec)
    atoms_amount == 0 && error("Standardizing cell failed!")

    Cell(clattice, cpositions, cnumbers)
end

find_primitive(cell::Cell; symprec::Real = 1e-5) = standardize_cell(cell; to_primitive = true, no_idealize = false, symprec = symprec)

refine_cell(cell::Cell; symprec::Real = 1e-5) = standardize_cell(cell; to_primitive = false, no_idealize = false, symprec = symprec)

function niggli_reduce(cell::Cell, symprec::Real = 1e-5)
    ccell = get_ccell(cell)
    clattice = getfield(ccell, :lattice)

    ret = ccall((:spg_niggli_reduce, spglib), Cint, (Ptr{Cdouble}, Cdouble), clattice, symprec)
    ret == 0 && error("Niggli reduce failed!")

    @set cell.lattice = clattice
end

function delaunay_reduce(cell::Cell, symprec::Real = 1e-5)
    ccell = get_ccell(cell)
    clattice = getfield(ccell, :lattice)

    ret = ccall((:spg_niggli_reduce, spglib), Cint, (Ptr{Cdouble}, Cdouble), clattice, symprec)
    ret == 0 && error("Delaunay reduce failed!")

    @set cell.lattice = clattice
end

function get_ir_reciprocal_mesh(cell::Cell,
                                grid::Vector{T},
                                shift::Vector{T} = [0, 0, 0],
                                is_time_reversal::Bool = true,
                                symprec::Real = 1e-5) where {T <: Integer}
    all(x -> x in (zero(T), one(T)), shift) || throw(ArgumentError("The shift can be only a vector of ones or zeros!"))

    qpoints_amount = prod(grid)
    grid_address = Array{Cint}(undef, qpoints_amount, 3)
    mapping = Array{Cint}(undef, qpoints_amount)
    ccell = get_ccell(cell)
    clattice, cpositions, cnumbers = getfields(ccell, :lattice, :positions, :numbers)

    ret = ccall((:spg_get_ir_reciprocal_mesh, spglib), Cint,
        (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        grid_address, mapping, grid, shift, is_time_reversal, clattice, cpositions, cnumbers, length(cnumbers), symprec)
    ret != qpoints_amount && error("Something wrong happens when finding mesh!")

    mapping, grid_address
end

function get_stabilized_reciprocal_mesh(rotations::Vector{Matrix{T}},
                                        grid::Vector{T},
                                        shift::Vector{T} = [0, 0, 0],
                                        qpoints::Vector{} = nothing,
                                        is_time_reversal::Bool = true) where {T <: Integer}
    all(x -> x in (zero(T), one(T)), shift) || throw(ArgumentError("The shift can be only a vector of ones or zeros!"))

    qpoints_amount = prod(grid)
    grid_address = Array{Cint}(undef, qpoints_amount, 3)
    mapping_table = Array{Cint}(undef, qpoints_amount)
    qpoints == nothing ? qpoints = Float64[0, 0, 0] : qpoints = Vector(qpoints)

    ret = ccall((:spg_get_stabilized_reciprocal_mesh, spglib), Cint,
        (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Ptr{Cint}, Cint, Ptr{Cint}),
        grid_address, mapping_table, grid, shift, is_time_reversal, length(rotations), rotations, length(qpoints), qpoints)
    ret != qpoints_amount && error("Something wrong happens when finding mesh!")

    mapping_table, grid_address
end

end
