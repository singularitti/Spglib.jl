# This is an internal type, do not export!
struct SpglibDataset
    spacegroup_number::Cint
    hall_number::Cint
    international_symbol::NTuple{11,Cchar}
    hall_symbol::NTuple{17,Cchar}
    choice::NTuple{6,Cchar}
    transformation_matrix::NTuple{9,Cdouble}
    origin_shift::NTuple{3,Cdouble}
    n_operations::Cint
    rotations::Ptr{NTuple{9,Cint}}
    translations::Ptr{NTuple{3,Cdouble}}
    n_atoms::Cint
    wyckoffs::Ptr{Cint}
    site_symmetry_symbols::Ptr{NTuple{7,Cchar}}
    equivalent_atoms::Ptr{Cint}
    mapping_to_primitive::Ptr{Cint}
    n_std_atoms::Cint
    std_lattice::NTuple{9,Cdouble}
    std_types::Ptr{Cint}
    std_positions::Ptr{NTuple{3,Cdouble}}
    std_rotation_matrix::NTuple{9,Cdouble}
    std_mapping_to_primitive::Ptr{Cint}
    pointgroup_symbol::NTuple{6,Cchar}
end

"""
    Dataset(spacegroup_number, hall_number, international_symbol, hall_symbol, choice, transformation_matrix, origin_shift, n_operations, rotations, translations, n_atoms, wyckoffs, site_symmetry_symbols, equivalent_atoms, mapping_to_primitive, n_std_atoms, std_lattice, std_types, std_positions, std_rotation_matrix, std_mapping_to_primitive, pointgroup_symbol)

Represent `SpglibDataset`, see its [official documentation](https://spglib.github.io/spglib/dataset.html#spglib-dataset).

!!! note
    Fields `crystallographic_orbits` and `primitive_lattice` are added after `spglib` `v1.15.0`.
"""
struct Dataset
    spacegroup_number::Int
    hall_number::Int
    international_symbol::String
    hall_symbol::String
    choice::String
    transformation_matrix::Matrix{Float64}
    origin_shift::Vector{Float64}
    n_operations::Int
    rotations::Array{Float64,3}
    translations::Matrix{Float64}
    n_atoms::Int
    wyckoffs::Vector{String}
    site_symmetry_symbols::Vector{String}
    equivalent_atoms::Vector{Int}
    mapping_to_primitive::Vector{Int}
    n_std_atoms::Int
    std_lattice::Matrix{Float64}
    std_types::Vector{Int}
    std_positions::Matrix{Float64}
    std_rotation_matrix::Matrix{Float64}
    std_mapping_to_primitive::Vector{Int}
    pointgroup_symbol::String
end

function Base.convert(::Type{Dataset}, dataset::SpglibDataset)
    r = unsafe_wrap(Vector{NTuple{9,Cint}}, dataset.rotations, dataset.n_operations)
    t = unsafe_wrap(Vector{NTuple{3,Float64}}, dataset.translations, dataset.n_operations)
    wyckoffs = unsafe_wrap(Vector{Cint}, dataset.wyckoffs, dataset.n_atoms)
    pos = unsafe_wrap(Vector{NTuple{3,Float64}}, dataset.std_positions, dataset.n_std_atoms)
    return Dataset(
        dataset.spacegroup_number,
        dataset.hall_number,
        cchars2string(dataset.international_symbol),
        cchars2string(dataset.hall_symbol),
        cchars2string(dataset.choice),
        tuple2matrix(dataset.transformation_matrix),
        collect(dataset.origin_shift),
        dataset.n_operations,
        rotsFromTuple(r, dataset.n_operations),
        transFromTuple(t, dataset.n_operations),
        dataset.n_atoms,
        [string(LETTERS[x+1]) for x in wyckoffs],  # Need to add 1 because of C-index starts from 0
        map(
            cchars2string,
            unsafe_wrap(
                Vector{NTuple{7,Cchar}},
                dataset.site_symmetry_symbols,
                dataset.n_atoms,
            ),
        ),
        unsafe_wrap(Vector{Cint}, dataset.equivalent_atoms, dataset.n_atoms),
        unsafe_wrap(Vector{Cint}, dataset.mapping_to_primitive, dataset.n_atoms),
        dataset.n_std_atoms,
        transpose(tuple2matrix(dataset.std_lattice)),
        unsafe_wrap(Vector{Cint}, dataset.std_types, dataset.n_std_atoms),
        transFromTuple(pos, dataset.n_std_atoms),
        tuple2matrix(dataset.std_rotation_matrix),
        unsafe_wrap(Vector{Cint}, dataset.std_mapping_to_primitive, dataset.n_std_atoms),
        cchars2string(dataset.pointgroup_symbol),
    )
end
