module Spglib

using spglib_jll: libsymspg

# All public methods
export get_symmetry,
    get_symmetry!,
    get_symmetry_with_collinear_spin,
    get_symmetry_with_collinear_spin!,
    get_symmetry_from_database,
    get_symmetry_from_database!,
    get_hall_number_from_symmetry,
    get_dataset,
    get_dataset_with_hall_number,
    get_spacegroup_number,
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
    get_stabilized_reciprocal_mesh,
    get_version

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

include("model.jl")
include("error.jl")
include("symmetry.jl")
include("standardize.jl")
include("reduce.jl")
include("reciprocal.jl")

end
