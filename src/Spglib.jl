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

include("model.jl")
include("symmetry.jl")
include("standardize.jl")
include("reduce.jl")
include("reciprocal.jl")

end
