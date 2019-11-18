# Automatically generated using Clang.jl


struct Cell
    size::Cint
    lattice::NTuple{3, NTuple{3, Cdouble}}
    types::Ptr{Cint}
    position::Ptr{Cvoid}
end

# Skipping MacroDefinition: debug_print ( ... )

# const debug_print_matrix_d3 = a
# const debug_print_matrix_i3 = a

# Skipping MacroDefinition: debug_print_vectors_d3 ( ... )
# Skipping MacroDefinition: debug_print_vectors_with_label ( ... )
# Skipping MacroDefinition: warning_print ( ... )

const KPT_NUM_BZ_SEARCH_SPACE = 125

struct MatINT
    size::Cint
    mat::Ptr{Cvoid}
end

struct VecDBL
    size::Cint
    vec::Ptr{Cvoid}
end

const NIGGLI_MAJOR_VERSION = 0
const NIGGLI_MINOR_VERSION = 1
const NIGGLI_MICRO_VERSION = 2

@cenum Holohedry::UInt32 begin
    HOLOHEDRY_NONE = 0
    TRICLI = 1
    MONOCLI = 2
    ORTHO = 3
    TETRA = 4
    TRIGO = 5
    HEXA = 6
    CUBIC = 7
end

@cenum Laue::UInt32 begin
    LAUE_NONE = 0
    LAUE1 = 1
    LAUE2M = 2
    LAUEMMM = 3
    LAUE4M = 4
    LAUE4MMM = 5
    LAUE3 = 6
    LAUE3M = 7
    LAUE6M = 8
    LAUE6MMM = 9
    LAUEM3 = 10
    LAUEM3M = 11
end


struct Pointgroup
    number::Cint
    symbol::NTuple{6, UInt8}
    schoenflies::NTuple{4, UInt8}
    holohedry::Holohedry
    laue::Laue
end

struct Primitive
    cell::Ptr{Cell}
    mapping_table::Ptr{Cint}
    size::Cint
    t_mat::NTuple{3, NTuple{3, Cdouble}}
    tolerance::Cdouble
end

struct Spacegroup
    number::Cint
    hall_number::Cint
    pointgroup_number::Cint
    schoenflies::NTuple{7, UInt8}
    hall_symbol::NTuple{17, UInt8}
    international::NTuple{32, UInt8}
    international_long::NTuple{20, UInt8}
    international_short::NTuple{11, UInt8}
    choice::NTuple{6, UInt8}
    bravais_lattice::NTuple{3, NTuple{3, Cdouble}}
    origin_shift::NTuple{3, Cdouble}
end

@cenum Centering::UInt32 begin
    CENTERING_ERROR = 0
    PRIMITIVE = 1
    BODY = 2
    FACE = 3
    A_FACE = 4
    B_FACE = 5
    C_FACE = 6
    BASE = 7
    R_CENTER = 8
end


struct SpacegroupType
    number::Cint
    schoenflies::NTuple{7, UInt8}
    hall_symbol::NTuple{17, UInt8}
    international::NTuple{32, UInt8}
    international_full::NTuple{20, UInt8}
    international_short::NTuple{11, UInt8}
    choice::NTuple{6, UInt8}
    centering::Centering
    pointgroup_number::Cint
end

struct SpglibDataset
    spacegroup_number::Cint
    hall_number::Cint
    international_symbol::NTuple{11, UInt8}
    hall_symbol::NTuple{17, UInt8}
    choice::NTuple{6, UInt8}
    transformation_matrix::NTuple{3, NTuple{3, Cdouble}}
    origin_shift::NTuple{3, Cdouble}
    n_operations::Cint
    rotations::Ptr{Cvoid}
    translations::Ptr{Cvoid}
    n_atoms::Cint
    wyckoffs::Ptr{Cint}
    equivalent_atoms::Ptr{Cint}
    n_std_atoms::Cint
    std_lattice::NTuple{3, NTuple{3, Cdouble}}
    std_types::Ptr{Cint}
    std_positions::Ptr{Cvoid}
    pointgroup_symbol::NTuple{6, UInt8}
end

struct SpglibSpacegroupType
    number::Cint
    international_short::NTuple{11, UInt8}
    international_full::NTuple{20, UInt8}
    international::NTuple{32, UInt8}
    schoenflies::NTuple{7, UInt8}
    hall_symbol::NTuple{17, UInt8}
    choice::NTuple{6, UInt8}
    pointgroup_international::NTuple{6, UInt8}
    pointgroup_schoenflies::NTuple{4, UInt8}
    arithmetic_crystal_class_number::Cint
    arithmetic_crystal_class_symbol::NTuple{7, UInt8}
end

struct Symmetry
    size::Cint
    rot::Ptr{Cvoid}
    trans::Ptr{Cvoid}
end

struct PointSymmetry
    rot::NTuple{48, NTuple{3, NTuple{3, Cint}}}
    size::Cint
end

const SPGLIB_MAJOR_VERSION = 1
const SPGLIB_MINOR_VERSION = 9
const SPGLIB_MICRO_VERSION = 4
