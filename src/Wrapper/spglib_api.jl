# Julia wrapper for header: arithmetic.h
# Automatically generated using Clang.jl


function arth_get_symbol(symbol, spgroup_number)
    ccall((:arth_get_symbol, arithmetic), Cint, (Ptr{UInt8}, Cint), symbol, spgroup_number)
end
# Julia wrapper for header: cell.h
# Automatically generated using Clang.jl


function cel_alloc_cell(size)
    ccall((:cel_alloc_cell, cell), Ptr{Cell}, (Cint,), size)
end

function cel_free_cell(cell)
    ccall((:cel_free_cell, cell), Cvoid, (Ptr{Cell},), cell)
end

function cel_set_cell(cell, lattice, position, types)
    ccall((:cel_set_cell, cell), Cvoid, (Ptr{Cell}, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}), cell, lattice, position, types)
end

function cel_copy_cell(cell)
    ccall((:cel_copy_cell, cell), Ptr{Cell}, (Ptr{Cell},), cell)
end

function cel_is_overlap(a, b, lattice, symprec)
    ccall((:cel_is_overlap, cell), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{NTuple{3, Cdouble}}, Cdouble), a, b, lattice, symprec)
end

function cel_trim_cell(mapping_table, trimmed_lattice, cell, symprec)
    ccall((:cel_trim_cell, cell), Ptr{Cell}, (Ptr{Cint}, Ptr{NTuple{3, Cdouble}}, Ptr{Cell}, Cdouble), mapping_table, trimmed_lattice, cell, symprec)
end
# Julia wrapper for header: debug.h
# Automatically generated using Clang.jl

# Julia wrapper for header: delaunay.h
# Automatically generated using Clang.jl


function del_delaunay_reduce(lattice_new, lattice, symprec)
    ccall((:del_delaunay_reduce, delaunay), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Cdouble), lattice_new, lattice, symprec)
end

function del_delaunay_reduce_2D(min_lattice, lattice, unique_axis, symprec)
    ccall((:del_delaunay_reduce_2D, delaunay), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Cint, Cdouble), min_lattice, lattice, unique_axis, symprec)
end
# Julia wrapper for header: hall_symbol.h
# Automatically generated using Clang.jl


function hal_match_hall_symbol_db(origin_shift, bravais_lattice, hall_number, centering, symmetry, symprec)
    ccall((:hal_match_hall_symbol_db, hall_symbol), Cint, (Ptr{Cdouble}, Ptr{NTuple{3, Cdouble}}, Cint, Centering, Ptr{Symmetry}, Cdouble), origin_shift, bravais_lattice, hall_number, centering, symmetry, symprec)
end
# Julia wrapper for header: kgrid.h
# Automatically generated using Clang.jl


function kgd_get_all_grid_addresses(grid_address, mesh)
    ccall((:kgd_get_all_grid_addresses, kgrid), Cvoid, (Ptr{NTuple{3, Cint}}, Ptr{Cint}), grid_address, mesh)
end

function kgd_get_grid_point_double_mesh(address_double, mesh)
    ccall((:kgd_get_grid_point_double_mesh, kgrid), Cint, (Ptr{Cint}, Ptr{Cint}), address_double, mesh)
end

function kgd_get_grid_address_double_mesh(address_double, address, mesh, is_shift)
    ccall((:kgd_get_grid_address_double_mesh, kgrid), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), address_double, address, mesh, is_shift)
end
# Julia wrapper for header: kpoint.h
# Automatically generated using Clang.jl


function kpt_get_irreducible_reciprocal_mesh(grid_address, map, mesh, is_shift, rot_reciprocal)
    ccall((:kpt_get_irreducible_reciprocal_mesh, kpoint), Cint, (Ptr{NTuple{3, Cint}}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{MatINT}), grid_address, map, mesh, is_shift, rot_reciprocal)
end

function kpt_get_stabilized_reciprocal_mesh(grid_address, map, mesh, is_shift, is_time_reversal, rotations, num_q, qpoints)
    ccall((:kpt_get_stabilized_reciprocal_mesh, kpoint), Cint, (Ptr{NTuple{3, Cint}}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Ptr{MatINT}, Cint, Ptr{NTuple{3, Cdouble}}), grid_address, map, mesh, is_shift, is_time_reversal, rotations, num_q, qpoints)
end

function kpt_get_grid_points_by_rotations(rot_grid_points, address_orig, rot_reciprocal, mesh, is_shift)
    ccall((:kpt_get_grid_points_by_rotations, kpoint), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{MatINT}, Ptr{Cint}, Ptr{Cint}), rot_grid_points, address_orig, rot_reciprocal, mesh, is_shift)
end

function kpt_get_BZ_grid_points_by_rotations(rot_grid_points, address_orig, rot_reciprocal, mesh, is_shift, bz_map)
    ccall((:kpt_get_BZ_grid_points_by_rotations, kpoint), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{MatINT}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), rot_grid_points, address_orig, rot_reciprocal, mesh, is_shift, bz_map)
end

function kpt_relocate_BZ_grid_address(bz_grid_address, bz_map, grid_address, mesh, rec_lattice, is_shift)
    ccall((:kpt_relocate_BZ_grid_address, kpoint), Cint, (Ptr{NTuple{3, Cint}}, Ptr{Cint}, Ptr{NTuple{3, Cint}}, Ptr{Cint}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}), bz_grid_address, bz_map, grid_address, mesh, rec_lattice, is_shift)
end

function kpt_get_point_group_reciprocal(rotations, is_time_reversal)
    ccall((:kpt_get_point_group_reciprocal, kpoint), Ptr{MatINT}, (Ptr{MatINT}, Cint), rotations, is_time_reversal)
end

function kpt_get_point_group_reciprocal_with_q(rot_reciprocal, symprec, num_q, qpoints)
    ccall((:kpt_get_point_group_reciprocal_with_q, kpoint), Ptr{MatINT}, (Ptr{MatINT}, Cdouble, Cint, Ptr{NTuple{3, Cdouble}}), rot_reciprocal, symprec, num_q, qpoints)
end
# Julia wrapper for header: mathfunc.h
# Automatically generated using Clang.jl


function mat_get_determinant_d3(a)
    ccall((:mat_get_determinant_d3, mathfunc), Cdouble, (Ptr{NTuple{3, Cdouble}},), a)
end

function mat_get_determinant_i3(a)
    ccall((:mat_get_determinant_i3, mathfunc), Cint, (Ptr{NTuple{3, Cint}},), a)
end

function mat_get_trace_i3(a)
    ccall((:mat_get_trace_i3, mathfunc), Cint, (Ptr{NTuple{3, Cint}},), a)
end

function mat_copy_matrix_d3(a, b)
    ccall((:mat_copy_matrix_d3, mathfunc), Cvoid, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}), a, b)
end

function mat_copy_matrix_i3(a, b)
    ccall((:mat_copy_matrix_i3, mathfunc), Cvoid, (Ptr{NTuple{3, Cint}}, Ptr{NTuple{3, Cint}}), a, b)
end

function mat_copy_vector_d3(a, b)
    ccall((:mat_copy_vector_d3, mathfunc), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}), a, b)
end

function mat_copy_vector_i3(a, b)
    ccall((:mat_copy_vector_i3, mathfunc), Cvoid, (Ptr{Cint}, Ptr{Cint}), a, b)
end

function mat_check_identity_matrix_i3(a, b)
    ccall((:mat_check_identity_matrix_i3, mathfunc), Cint, (Ptr{NTuple{3, Cint}}, Ptr{NTuple{3, Cint}}), a, b)
end

function mat_check_identity_matrix_d3(a, b, symprec)
    ccall((:mat_check_identity_matrix_d3, mathfunc), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Cdouble), a, b, symprec)
end

function mat_check_identity_matrix_id3(a, b, symprec)
    ccall((:mat_check_identity_matrix_id3, mathfunc), Cint, (Ptr{NTuple{3, Cint}}, Ptr{NTuple{3, Cdouble}}, Cdouble), a, b, symprec)
end

function mat_multiply_matrix_d3(m, a, b)
    ccall((:mat_multiply_matrix_d3, mathfunc), Cvoid, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}), m, a, b)
end

function mat_multiply_matrix_i3(m, a, b)
    ccall((:mat_multiply_matrix_i3, mathfunc), Cvoid, (Ptr{NTuple{3, Cint}}, Ptr{NTuple{3, Cint}}, Ptr{NTuple{3, Cint}}), m, a, b)
end

function mat_multiply_matrix_di3(m, a, b)
    ccall((:mat_multiply_matrix_di3, mathfunc), Cvoid, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cint}}), m, a, b)
end

function mat_multiply_matrix_id3(m, a, b)
    ccall((:mat_multiply_matrix_id3, mathfunc), Cvoid, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cint}}, Ptr{NTuple{3, Cdouble}}), m, a, b)
end

function mat_multiply_matrix_vector_i3(v, a, b)
    ccall((:mat_multiply_matrix_vector_i3, mathfunc), Cvoid, (Ptr{Cint}, Ptr{NTuple{3, Cint}}, Ptr{Cint}), v, a, b)
end

function mat_multiply_matrix_vector_d3(v, a, b)
    ccall((:mat_multiply_matrix_vector_d3, mathfunc), Cvoid, (Ptr{Cdouble}, Ptr{NTuple{3, Cdouble}}, Ptr{Cdouble}), v, a, b)
end

function mat_multiply_matrix_vector_id3(v, a, b)
    ccall((:mat_multiply_matrix_vector_id3, mathfunc), Cvoid, (Ptr{Cdouble}, Ptr{NTuple{3, Cint}}, Ptr{Cdouble}), v, a, b)
end

function mat_multiply_matrix_vector_di3(v, a, b)
    ccall((:mat_multiply_matrix_vector_di3, mathfunc), Cvoid, (Ptr{Cdouble}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}), v, a, b)
end

function mat_add_matrix_i3(m, a, b)
    ccall((:mat_add_matrix_i3, mathfunc), Cvoid, (Ptr{NTuple{3, Cint}}, Ptr{NTuple{3, Cint}}, Ptr{NTuple{3, Cint}}), m, a, b)
end

function mat_cast_matrix_3i_to_3d(m, a)
    ccall((:mat_cast_matrix_3i_to_3d, mathfunc), Cvoid, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cint}}), m, a)
end

function mat_cast_matrix_3d_to_3i(m, a)
    ccall((:mat_cast_matrix_3d_to_3i, mathfunc), Cvoid, (Ptr{NTuple{3, Cint}}, Ptr{NTuple{3, Cdouble}}), m, a)
end

function mat_inverse_matrix_d3(m, a, precision)
    ccall((:mat_inverse_matrix_d3, mathfunc), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Cdouble), m, a, precision)
end

function mat_get_similar_matrix_d3(m, a, b, precision)
    ccall((:mat_get_similar_matrix_d3, mathfunc), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Cdouble), m, a, b, precision)
end

function mat_transpose_matrix_d3(a, b)
    ccall((:mat_transpose_matrix_d3, mathfunc), Cvoid, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}), a, b)
end

function mat_transpose_matrix_i3(a, b)
    ccall((:mat_transpose_matrix_i3, mathfunc), Cvoid, (Ptr{NTuple{3, Cint}}, Ptr{NTuple{3, Cint}}), a, b)
end

function mat_get_metric(metric, lattice)
    ccall((:mat_get_metric, mathfunc), Cvoid, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}), metric, lattice)
end

function mat_norm_squared_d3(a)
    ccall((:mat_norm_squared_d3, mathfunc), Cdouble, (Ptr{Cdouble},), a)
end

function mat_norm_squared_i3(a)
    ccall((:mat_norm_squared_i3, mathfunc), Cint, (Ptr{Cint},), a)
end

function mat_Dabs(a)
    ccall((:mat_Dabs, mathfunc), Cdouble, (Cdouble,), a)
end

function mat_Nint(a)
    ccall((:mat_Nint, mathfunc), Cint, (Cdouble,), a)
end

function mat_Dmod1(a)
    ccall((:mat_Dmod1, mathfunc), Cdouble, (Cdouble,), a)
end

function mat_alloc_MatINT(size)
    ccall((:mat_alloc_MatINT, mathfunc), Ptr{MatINT}, (Cint,), size)
end

function mat_free_MatINT(matint)
    ccall((:mat_free_MatINT, mathfunc), Cvoid, (Ptr{MatINT},), matint)
end

function mat_alloc_VecDBL(size)
    ccall((:mat_alloc_VecDBL, mathfunc), Ptr{VecDBL}, (Cint,), size)
end

function mat_free_VecDBL(vecdbl)
    ccall((:mat_free_VecDBL, mathfunc), Cvoid, (Ptr{VecDBL},), vecdbl)
end

function mat_is_int_matrix(mat, symprec)
    ccall((:mat_is_int_matrix, mathfunc), Cint, (Ptr{NTuple{3, Cdouble}}, Cdouble), mat, symprec)
end
# Julia wrapper for header: niggli.h
# Automatically generated using Clang.jl


function niggli_get_major_version()
    ccall((:niggli_get_major_version, niggli), Cint, ())
end

function niggli_get_minor_version()
    ccall((:niggli_get_minor_version, niggli), Cint, ())
end

function niggli_get_micro_version()
    ccall((:niggli_get_micro_version, niggli), Cint, ())
end

function niggli_reduce(lattice_, eps_)
    ccall((:niggli_reduce, niggli), Cint, (Ptr{Cdouble}, Cdouble), lattice_, eps_)
end
# Julia wrapper for header: pointgroup.h
# Automatically generated using Clang.jl


function ptg_get_transformation_matrix(transform_mat, rotations, num_rotations)
    ccall((:ptg_get_transformation_matrix, pointgroup), Pointgroup, (Ptr{NTuple{3, Cint}}, Ptr{NTuple{3, NTuple{3, Cint}}}, Cint), transform_mat, rotations, num_rotations)
end

function ptg_get_pointgroup(pointgroup_number)
    ccall((:ptg_get_pointgroup, pointgroup), Pointgroup, (Cint,), pointgroup_number)
end

function ptg_get_pointsymmetry(rotations, num_rotations)
    ccall((:ptg_get_pointsymmetry, pointgroup), PointSymmetry, (Ptr{NTuple{3, NTuple{3, Cint}}}, Cint), rotations, num_rotations)
end
# Julia wrapper for header: primitive.h
# Automatically generated using Clang.jl


function prm_alloc_primitive(size)
    ccall((:prm_alloc_primitive, primitive), Ptr{Primitive}, (Cint,), size)
end

function prm_free_primitive(primitive)
    ccall((:prm_free_primitive, primitive), Cvoid, (Ptr{Primitive},), primitive)
end

function prm_get_primitive(cell, symprec)
    ccall((:prm_get_primitive, primitive), Ptr{Primitive}, (Ptr{Cell}, Cdouble), cell, symprec)
end
# Julia wrapper for header: refinement.h
# Automatically generated using Clang.jl


function ref_get_refined_symmetry_operations(cell, primitive, spacegroup, symprec)
    ccall((:ref_get_refined_symmetry_operations, refinement), Ptr{Symmetry}, (Ptr{Cell}, Ptr{Cell}, Ptr{Spacegroup}, Cdouble), cell, primitive, spacegroup, symprec)
end

function ref_get_Wyckoff_positions(wyckoffs, equiv_atoms, primitive, cell, spacegroup, symmetry, mapping_table, symprec)
    ccall((:ref_get_Wyckoff_positions, refinement), Ptr{Cell}, (Ptr{Cint}, Ptr{Cint}, Ptr{Cell}, Ptr{Cell}, Ptr{Spacegroup}, Ptr{Symmetry}, Ptr{Cint}, Cdouble), wyckoffs, equiv_atoms, primitive, cell, spacegroup, symmetry, mapping_table, symprec)
end
# Julia wrapper for header: site_symmetry.h
# Automatically generated using Clang.jl


function ssm_get_exact_positions(wyckoffs, equiv_atoms, bravais, conv_sym, hall_number, symprec)
    ccall((:ssm_get_exact_positions, site_symmetry), Ptr{VecDBL}, (Ptr{Cint}, Ptr{Cint}, Ptr{Cell}, Ptr{Symmetry}, Cint, Cdouble), wyckoffs, equiv_atoms, bravais, conv_sym, hall_number, symprec)
end
# Julia wrapper for header: sitesym_database.h
# Automatically generated using Clang.jl


function ssmdb_get_coordinate(rot, trans, index)
    ccall((:ssmdb_get_coordinate, sitesym_database), Cint, (Ptr{NTuple{3, Cint}}, Ptr{Cdouble}, Cint), rot, trans, index)
end

function ssmdb_get_wyckoff_indices(indices, index)
    ccall((:ssmdb_get_wyckoff_indices, sitesym_database), Cvoid, (Ptr{Cint}, Cint), indices, index)
end
# Julia wrapper for header: spacegroup.h
# Automatically generated using Clang.jl


function spa_get_spacegroup(spacegroup, cell, symprec)
    ccall((:spa_get_spacegroup, spacegroup), Ptr{Primitive}, (Ptr{Spacegroup}, Ptr{Cell}, Cdouble), spacegroup, cell, symprec)
end

function spa_get_spacegroup_with_hall_number(primitive, hall_number)
    ccall((:spa_get_spacegroup_with_hall_number, spacegroup), Spacegroup, (Ptr{Primitive}, Cint), primitive, hall_number)
end

function spa_transform_to_primitive(cell, trans_mat, centering, symprec)
    ccall((:spa_transform_to_primitive, spacegroup), Ptr{Cell}, (Ptr{Cell}, Ptr{NTuple{3, Cdouble}}, Centering, Cdouble), cell, trans_mat, centering, symprec)
end
# Julia wrapper for header: spg_database.h
# Automatically generated using Clang.jl


function spgdb_get_operation(rot, trans, hall_number)
    ccall((:spgdb_get_operation, spg_database), Cint, (Ptr{NTuple{3, Cint}}, Ptr{Cdouble}, Cint), rot, trans, hall_number)
end

function spgdb_get_operation_index(indices, hall_number)
    ccall((:spgdb_get_operation_index, spg_database), Cvoid, (Ptr{Cint}, Cint), indices, hall_number)
end

function spgdb_get_spacegroup_operations(hall_number)
    ccall((:spgdb_get_spacegroup_operations, spg_database), Ptr{Symmetry}, (Cint,), hall_number)
end

function spgdb_get_spacegroup_type(hall_number)
    ccall((:spgdb_get_spacegroup_type, spg_database), SpacegroupType, (Cint,), hall_number)
end
# Julia wrapper for header: spglib.h
# Automatically generated using Clang.jl


function spg_get_major_version()
    ccall((:spg_get_major_version, spglib), Cint, ())
end

function spg_get_minor_version()
    ccall((:spg_get_minor_version, spglib), Cint, ())
end

function spg_get_micro_version()
    ccall((:spg_get_micro_version, spglib), Cint, ())
end

function spg_get_dataset(lattice, position, types, num_atom, symprec)
    ccall((:spg_get_dataset, spglib), Ptr{SpglibDataset}, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble), lattice, position, types, num_atom, symprec)
end

function spgat_get_dataset(lattice, position, types, num_atom, symprec, angle_tolerance)
    ccall((:spgat_get_dataset, spglib), Ptr{SpglibDataset}, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble, Cdouble), lattice, position, types, num_atom, symprec, angle_tolerance)
end

function spg_get_dataset_with_hall_number(lattice, position, types, num_atom, hall_number, symprec)
    ccall((:spg_get_dataset_with_hall_number, spglib), Ptr{SpglibDataset}, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cint, Cdouble), lattice, position, types, num_atom, hall_number, symprec)
end

function spgat_get_dataset_with_hall_number(lattice, position, types, num_atom, hall_number, symprec, angle_tolerance)
    ccall((:spgat_get_dataset_with_hall_number, spglib), Ptr{SpglibDataset}, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cint, Cdouble, Cdouble), lattice, position, types, num_atom, hall_number, symprec, angle_tolerance)
end

function spg_free_dataset(dataset)
    ccall((:spg_free_dataset, spglib), Cvoid, (Ptr{SpglibDataset},), dataset)
end

function spg_get_symmetry(rotation, translation, max_size, lattice, position, types, num_atom, symprec)
    ccall((:spg_get_symmetry, spglib), Cint, (Ptr{NTuple{3, NTuple{3, Cint}}}, Ptr{NTuple{3, Cdouble}}, Cint, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble), rotation, translation, max_size, lattice, position, types, num_atom, symprec)
end

function spgat_get_symmetry(rotation, translation, max_size, lattice, position, types, num_atom, symprec, angle_tolerance)
    ccall((:spgat_get_symmetry, spglib), Cint, (Ptr{NTuple{3, NTuple{3, Cint}}}, Ptr{NTuple{3, Cdouble}}, Cint, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble, Cdouble), rotation, translation, max_size, lattice, position, types, num_atom, symprec, angle_tolerance)
end

function spg_get_symmetry_numerical(rotation, translation, max_size, lattice, position, types, num_atom, symprec)
    ccall((:spg_get_symmetry_numerical, spglib), Cint, (Ptr{NTuple{3, NTuple{3, Cint}}}, Ptr{NTuple{3, Cdouble}}, Cint, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble), rotation, translation, max_size, lattice, position, types, num_atom, symprec)
end

function spgat_get_symmetry_numerical(rotation, translation, max_size, lattice, position, types, num_atom, symprec, angle_tolerance)
    ccall((:spgat_get_symmetry_numerical, spglib), Cint, (Ptr{NTuple{3, NTuple{3, Cint}}}, Ptr{NTuple{3, Cdouble}}, Cint, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble, Cdouble), rotation, translation, max_size, lattice, position, types, num_atom, symprec, angle_tolerance)
end

function spg_get_symmetry_with_collinear_spin(rotation, translation, equivalent_atoms, max_size, lattice, position, types, spins, num_atom, symprec)
    ccall((:spg_get_symmetry_with_collinear_spin, spglib), Cint, (Ptr{NTuple{3, NTuple{3, Cint}}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Ptr{Cdouble}, Cint, Cdouble), rotation, translation, equivalent_atoms, max_size, lattice, position, types, spins, num_atom, symprec)
end

function spgat_get_symmetry_with_collinear_spin(rotation, translation, equivalent_atoms, max_size, lattice, position, types, spins, num_atom, symprec, angle_tolerance)
    ccall((:spgat_get_symmetry_with_collinear_spin, spglib), Cint, (Ptr{NTuple{3, NTuple{3, Cint}}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Ptr{Cdouble}, Cint, Cdouble, Cdouble), rotation, translation, equivalent_atoms, max_size, lattice, position, types, spins, num_atom, symprec, angle_tolerance)
end

function spg_get_multiplicity(lattice, position, types, num_atom, symprec)
    ccall((:spg_get_multiplicity, spglib), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble), lattice, position, types, num_atom, symprec)
end

function spgat_get_multiplicity(lattice, position, types, num_atom, symprec, angle_tolerance)
    ccall((:spgat_get_multiplicity, spglib), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble, Cdouble), lattice, position, types, num_atom, symprec, angle_tolerance)
end

function spg_get_international(symbol, lattice, position, types, num_atom, symprec)
    ccall((:spg_get_international, spglib), Cint, (Ptr{UInt8}, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble), symbol, lattice, position, types, num_atom, symprec)
end

function spgat_get_international(symbol, lattice, position, types, num_atom, symprec, angle_tolerance)
    ccall((:spgat_get_international, spglib), Cint, (Ptr{UInt8}, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble, Cdouble), symbol, lattice, position, types, num_atom, symprec, angle_tolerance)
end

function spg_get_schoenflies(symbol, lattice, position, types, num_atom, symprec)
    ccall((:spg_get_schoenflies, spglib), Cint, (Ptr{UInt8}, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble), symbol, lattice, position, types, num_atom, symprec)
end

function spgat_get_schoenflies(symbol, lattice, position, types, num_atom, symprec, angle_tolerance)
    ccall((:spgat_get_schoenflies, spglib), Cint, (Ptr{UInt8}, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble, Cdouble), symbol, lattice, position, types, num_atom, symprec, angle_tolerance)
end

function spg_get_pointgroup(symbol, trans_mat, rotations, num_rotations)
    ccall((:spg_get_pointgroup, spglib), Cint, (Ptr{UInt8}, Ptr{NTuple{3, Cint}}, Ptr{NTuple{3, NTuple{3, Cint}}}, Cint), symbol, trans_mat, rotations, num_rotations)
end

function spg_get_symmetry_from_database(rotations, translations, hall_number)
    ccall((:spg_get_symmetry_from_database, spglib), Cint, (Ptr{NTuple{3, NTuple{3, Cint}}}, Ptr{NTuple{3, Cdouble}}, Cint), rotations, translations, hall_number)
end

function spg_get_spacegroup_type(hall_number)
    ccall((:spg_get_spacegroup_type, spglib), SpglibSpacegroupType, (Cint,), hall_number)
end

function spg_standardize_cell(lattice, position, types, num_atom, to_primitive, no_idealize, symprec)
    ccall((:spg_standardize_cell, spglib), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cint, Cint, Cdouble), lattice, position, types, num_atom, to_primitive, no_idealize, symprec)
end

function spgat_standardize_cell(lattice, position, types, num_atom, to_primitive, no_idealize, symprec, angle_tolerance)
    ccall((:spgat_standardize_cell, spglib), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cint, Cint, Cdouble, Cdouble), lattice, position, types, num_atom, to_primitive, no_idealize, symprec, angle_tolerance)
end

function spg_find_primitive(lattice, position, types, num_atom, symprec)
    ccall((:spg_find_primitive, spglib), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble), lattice, position, types, num_atom, symprec)
end

function spgat_find_primitive(lattice, position, types, num_atom, symprec, angle_tolerance)
    ccall((:spgat_find_primitive, spglib), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble, Cdouble), lattice, position, types, num_atom, symprec, angle_tolerance)
end

function spg_refine_cell(lattice, position, types, num_atom, symprec)
    ccall((:spg_refine_cell, spglib), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble), lattice, position, types, num_atom, symprec)
end

function spgat_refine_cell(lattice, position, types, num_atom, symprec, angle_tolerance)
    ccall((:spgat_refine_cell, spglib), Cint, (Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble, Cdouble), lattice, position, types, num_atom, symprec, angle_tolerance)
end

function spg_delaunay_reduce(lattice, symprec)
    ccall((:spg_delaunay_reduce, spglib), Cint, (Ptr{NTuple{3, Cdouble}}, Cdouble), lattice, symprec)
end

function spg_get_grid_point_from_address(grid_address, mesh)
    ccall((:spg_get_grid_point_from_address, spglib), Cint, (Ptr{Cint}, Ptr{Cint}), grid_address, mesh)
end

function spg_get_ir_reciprocal_mesh(grid_address, map, mesh, is_shift, is_time_reversal, lattice, position, types, num_atom, symprec)
    ccall((:spg_get_ir_reciprocal_mesh, spglib), Cint, (Ptr{NTuple{3, Cint}}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Ptr{NTuple{3, Cdouble}}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}, Cint, Cdouble), grid_address, map, mesh, is_shift, is_time_reversal, lattice, position, types, num_atom, symprec)
end

function spg_get_stabilized_reciprocal_mesh(grid_address, map, mesh, is_shift, is_time_reversal, num_rot, rotations, num_q, qpoints)
    ccall((:spg_get_stabilized_reciprocal_mesh, spglib), Cint, (Ptr{NTuple{3, Cint}}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Ptr{NTuple{3, NTuple{3, Cint}}}, Cint, Ptr{NTuple{3, Cdouble}}), grid_address, map, mesh, is_shift, is_time_reversal, num_rot, rotations, num_q, qpoints)
end

function spg_get_grid_points_by_rotations(rot_grid_points, address_orig, num_rot, rot_reciprocal, mesh, is_shift)
    ccall((:spg_get_grid_points_by_rotations, spglib), Cint, (Ptr{Cint}, Ptr{Cint}, Cint, Ptr{NTuple{3, NTuple{3, Cint}}}, Ptr{Cint}, Ptr{Cint}), rot_grid_points, address_orig, num_rot, rot_reciprocal, mesh, is_shift)
end

function spg_get_BZ_grid_points_by_rotations(rot_grid_points, address_orig, num_rot, rot_reciprocal, mesh, is_shift, bz_map)
    ccall((:spg_get_BZ_grid_points_by_rotations, spglib), Cint, (Ptr{Cint}, Ptr{Cint}, Cint, Ptr{NTuple{3, NTuple{3, Cint}}}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), rot_grid_points, address_orig, num_rot, rot_reciprocal, mesh, is_shift, bz_map)
end

function spg_relocate_BZ_grid_address(bz_grid_address, bz_map, grid_address, mesh, rec_lattice, is_shift)
    ccall((:spg_relocate_BZ_grid_address, spglib), Cint, (Ptr{NTuple{3, Cint}}, Ptr{Cint}, Ptr{NTuple{3, Cint}}, Ptr{Cint}, Ptr{NTuple{3, Cdouble}}, Ptr{Cint}), bz_grid_address, bz_map, grid_address, mesh, rec_lattice, is_shift)
end

function spg_get_neighboring_grid_points(relative_grid_points, grid_point, relative_grid_address, num_relative_grid_address, mesh, bz_grid_address, bz_map)
    ccall((:spg_get_neighboring_grid_points, spglib), Cvoid, (Ptr{Cint}, Cint, Ptr{NTuple{3, Cint}}, Cint, Ptr{Cint}, Ptr{NTuple{3, Cint}}, Ptr{Cint}), relative_grid_points, grid_point, relative_grid_address, num_relative_grid_address, mesh, bz_grid_address, bz_map)
end

function spg_niggli_reduce(lattice, symprec)
    ccall((:spg_niggli_reduce, spglib), Cint, (Ptr{NTuple{3, Cdouble}}, Cdouble), lattice, symprec)
end
# Julia wrapper for header: spin.h
# Automatically generated using Clang.jl


function spn_get_collinear_operations(equiv_atoms, sym_nonspin, cell, spins, symprec)
    ccall((:spn_get_collinear_operations, spin), Ptr{Symmetry}, (Ptr{Cint}, Ptr{Symmetry}, Ptr{Cell}, Ptr{Cdouble}, Cdouble), equiv_atoms, sym_nonspin, cell, spins, symprec)
end
# Julia wrapper for header: symmetry.h
# Automatically generated using Clang.jl


function sym_alloc_symmetry(size)
    ccall((:sym_alloc_symmetry, symmetry), Ptr{Symmetry}, (Cint,), size)
end

function sym_free_symmetry(symmetry)
    ccall((:sym_free_symmetry, symmetry), Cvoid, (Ptr{Symmetry},), symmetry)
end

function sym_get_operation(primitive, symprec)
    ccall((:sym_get_operation, symmetry), Ptr{Symmetry}, (Ptr{Cell}, Cdouble), primitive, symprec)
end

function sym_reduce_operation(primitive, symmetry, symprec)
    ccall((:sym_reduce_operation, symmetry), Ptr{Symmetry}, (Ptr{Cell}, Ptr{Symmetry}, Cdouble), primitive, symmetry, symprec)
end

function sym_get_pure_translation(cell, symprec)
    ccall((:sym_get_pure_translation, symmetry), Ptr{VecDBL}, (Ptr{Cell}, Cdouble), cell, symprec)
end

function sym_reduce_pure_translation(cell, pure_trans, symprec)
    ccall((:sym_reduce_pure_translation, symmetry), Ptr{VecDBL}, (Ptr{Cell}, Ptr{VecDBL}, Cdouble), cell, pure_trans, symprec)
end

function sym_set_angle_tolerance(tolerance)
    ccall((:sym_set_angle_tolerance, symmetry), Cvoid, (Cdouble,), tolerance)
end

function sym_get_angle_tolerance()
    ccall((:sym_get_angle_tolerance, symmetry), Cdouble, ())
end
# Julia wrapper for header: version.h
# Automatically generated using Clang.jl

