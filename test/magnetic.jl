# See https://github.com/singularitti/Spglib.jl/issues/91#issuecomment-1206106977
@testset "Test example given by Jae-Mo Lihm (@jaemolihm)" begin
    lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    positions = [[-0.1, -0.1, -0.1], [0.1, 0.1, 0.1]]
    atoms = [1, 1]
    magmoms = [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]
    cell = Cell(lattice, positions, atoms, magmoms)
    rotations, translations = get_symmetry(cell, 1e-5)
    @test size(rotations) == size(translations) == (12,)
    rotations, translations, equivalent_atoms = get_magnetic_symmetry(cell, 1e-5)
    @test size(rotations) == size(translations) == (4,)
end

# From https://github.com/unkcpz/LibSymspg.jl/blob/53d2f6d/test/test_api.jl#L34-L77
@testset "Get symmetry operations" begin
    @testset "Normal symmetry" begin
        lattice = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
        positions = [[0, 0, 0], [0.5, 0.5, 0.5]]
        atoms = [1, 1]
        cell = Cell(lattice, positions, atoms, [0, 0])
        rotations, translations = get_symmetry(cell, 1e-5)
        @test size(rotations) == (96,)
        @test size(translations) == (96,)
        @test get_hall_number_from_symmetry(cell, 1e-5) == 529
    end
    # See https://github.com/spglib/spglib/blob/378240e/python/test/test_collinear_spin.py#L18-L37
    @testset "Get symmetry with collinear spins" begin
        lattice = [
            4.0 0.0 0.0
            0.0 4.0 0.0
            0.0 0.0 4.0
        ]
        positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        atoms = [1, 1]
        @testset "Test ferromagnetism" begin
            magmoms = [1.0, 1.0]
            cell = Cell(lattice, positions, atoms, magmoms)
            rotations, translations, equivalent_atoms = get_symmetry_with_collinear_spin(
                cell, 1e-5
            )
            @test size(rotations) == (96,)
            @test size(translations) == (96,)
            @test all(iszero(translation) for translation in translations[1:48])
            @test all(
                translation == [1 / 2, 1 / 2, 1 / 2] for translation in translations[49:96]
            )  # Compared with Python
            @test equivalent_atoms == [0, 0]
        end
        @testset "Test antiferromagnetism" begin
            magmoms = [1.0, -1.0]
            cell = Cell(lattice, positions, atoms, magmoms)
            rotations, translations, equivalent_atoms = get_symmetry_with_collinear_spin(
                cell, 1e-5
            )
            @test size(rotations) == (3, 3, 96)
            @test equivalent_atoms == [0, 0]
        end
        @testset "Test broken magmoms" begin
            magmoms = [1.0, 2.0]
            cell = Cell(lattice, positions, atoms, magmoms)
            rotations, translations, equivalent_atoms = get_symmetry_with_collinear_spin(
                cell, 1e-5
            )
            @test size(rotations) == (3, 3, 48)
            @test size(translations) == (3, 48)
            @test equivalent_atoms == [0, 1]
        end
    end
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_magnetic_dataset.py#L9-L44
@testset "Test Type-I" begin
    lattice = [
        6.8083 0.0 0.0
        0.0 6.8083 0.0
        0.0 0.0 12.3795
    ]
    positions = [
        [0.87664, 0.35295, 0.13499],
        [0.14705, 0.37664, 0.38499],
        [0.85295, 0.62336, 0.88499],
        [0.37664, 0.14705, 0.61501],
        [0.62336, 0.85295, 0.11501],
        [0.12336, 0.64705, 0.63499],
        [0.35295, 0.87664, 0.86501],
        [0.64705, 0.12336, 0.36501],
    ]
    atoms = [0, 0, 0, 0, 0, 0, 0, 0]
    magmoms = [
        [1.67, -8.9, 0.0],
        [8.9, 1.67, 0.0],
        [-8.9, -1.67, 0.0],
        [1.67, 8.9, 0.0],
        [-1.67, -8.9, 0.0],
        [-1.67, 8.9, 0.0],
        [-8.9, 1.67, 0.0],
        [8.9, -1.67, 0.0],
    ]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell, 1e-5)
    @test dataset.uni_number == 771
    @test dataset.msg_type == 1
    @test dataset.hall_number == 369
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 8
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [0 -1 0; 1 0 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [0 1 0; -1 0 0; 0 0 1],
        [1 0 0; 0 -1 0; 0 0 -1],
        [0 -1 0; -1 0 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [0 1 0; 1 0 0; 0 0 -1],
    ]  # Compared with Python results
    @test dataset.translations ≈ [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.25],
        [-1.11022302e-16, 1.11022302e-16, 0.5],
        [0.5, 0.5, 0.75],
        [0.5, 0.5, 0.75],
        [0.0, 0.0, 0.5],
        [0.5, 0.5, 0.25],
        [-1.11022302e-16, 1.11022302e-16, 0.0],
    ]  # Compared with Python results
    @test dataset.time_reversals == falses(8)
    @test dataset.n_atoms == 8
    @test dataset.equivalent_atoms == zeros(Int32, 8) .+ 1
    @test dataset.transformation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.origin_shift == [5.55111512e-17, -1.66533454e-16, 0]
    @test dataset.n_std_atoms == 8
    @test dataset.std_lattice ≈ Lattice([
        6.8083 0 0
        0 6.8083 0
        0 0 12.3795
    ])  # Compared with Python results
    @test dataset.std_types == fill(1, 8)
    @test dataset.std_positions ≈ [
        [0.87664, 0.35295, 0.13499],
        [0.14705, 0.37664, 0.38499],
        [0.85295, 0.62336, 0.88499],
        [0.37664, 0.14705, 0.61501],
        [0.62336, 0.85295, 0.11501],
        [0.12336, 0.64705, 0.63499],
        [0.35295, 0.87664, 0.86501],
        [0.64705, 0.12336, 0.36501],
    ]
    @test dataset.std_tensors == [
        [1.67, -8.9, 0.0],
        [8.9, 1.67, 0.0],
        [-8.9, -1.67, 0.0],
        [1.67, 8.9, 0.0],
        [-1.67, -8.9, 0.0],
        [-1.67, 8.9, 0.0],
        [-8.9, 1.67, 0.0],
        [8.9, -1.67, 0.0],
    ]
    @test dataset.std_rotation_matrix ≈ [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.primitive_lattice == Lattice([
        6.8083 0 0
        0 6.8083 0
        0 0 12.3795
    ])
end

# From https://github.com/spglib/spglib/blob/f6abb97/test/functional/fortran/test_fortran_spg_get_symmetry_with_site_tensors.F90#L46-L97
@testset "Test site tensors for rutile (type III)" begin
    lattice = [
        4.0 0.0 0.0
        0.0 4.0 0.0
        0.0 0.0 3.0
    ]
    positions =
        [
            0.0 0.0 0.0
            0.5 0.5 0.5
            0.3 0.3 0.0
            0.7 0.7 0.0
            0.2 0.8 0.5
            0.8 0.2 0.5
        ] .+ [0.1 0.1 0.0]
    atoms = [1, 1, 2, 2, 2, 2]
    magmoms = [1.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 1158
    @test dataset.msg_type == 3
    @test dataset.hall_number == 419
    @test dataset.tensor_rank == 0
    @test dataset.n_operations == 16
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [0 -1 0; 1 0 0; 0 0 1],
        [0 1 0; -1 0 0; 0 0 -1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
        [0 1 0; -1 0 0; 0 0 1],
        [0 -1 0; 1 0 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 1],
        [0 -1 0; -1 0 0; 0 0 -1],
        [0 1 0; 1 0 0; 0 0 1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [0 1 0; 1 0 0; 0 0 -1],
        [0 -1 0; -1 0 0; 0 0 1],
    ]  # Compared with Python results
    @test dataset.translations ≈ [
        [0.0, 0.0, 0.0],
        [0.2, 0.2, 0.0],
        [0.7, 0.5, 0.5],
        [0.5, 0.7, 0.5],
        [0.2, 0.2, 0.0],
        [0.0, 0.0, 0.0],
        [0.5, 0.7, 0.5],
        [0.7, 0.5, 0.5],
        [0.5, 0.7, 0.5],
        [0.7, 0.5, 0.5],
        [0.2, 0.2, 0.0],
        [0.0, 0.0, 0.0],
        [0.7, 0.5, 0.5],
        [0.5, 0.7, 0.5],
        [0.0, 0.0, 0.0],
        [0.2, 0.2, 0.0],
    ]  # Compared with Python results
    @test dataset.time_reversals == [
        false,
        false,
        true,
        true,
        false,
        false,
        true,
        true,
        true,
        true,
        false,
        false,
        true,
        true,
        false,
        false,
    ]  # Compared with Python results
    @test dataset.n_atoms == 6
    @test dataset.equivalent_atoms == [0, 0, 2, 2, 2, 2] .+ 1  # Compared with Python results
    @test dataset.transformation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.origin_shift == [0.9, 0.9, 0.0]
    @test dataset.n_std_atoms == 6
    @test dataset.std_lattice == Lattice([4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 3.0])
    @test dataset.std_types == [1, 1, 2, 2, 2, 2]
    @test dataset.std_positions ≈ [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.3, 0.3, 0.0],
        [0.7, 0.7, 0.0],
        [0.2, 0.8, 0.5],
        [0.8, 0.2, 0.5],
    ]
    @test dataset.std_tensors == [1.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    @test dataset.std_rotation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.primitive_lattice ==
        Lattice([0.0, 0.0, 3.0], [4.0, 0.0, 0.0], [0.0, 4.0, 0.0])
end

# From https://github.com/spglib/spglib/blob/f6abb97/test/functional/fortran/test_fortran_spg_get_symmetry_with_site_tensors.F90#L99-L146
@testset "Test site tensors for Cr (type IV)" begin
    lattice = Lattice([
        4.0 0.0 0.0
        0.0 4.0 0.0
        0.0 0.0 4.0
    ])
    positions = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5]] .+ Ref([0.1, 0.1, 0.0])
    atoms = [1, 1]
    magmoms = [1.0, -1.0]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 1009
    @test dataset.msg_type == 4
    @test dataset.hall_number == 400
    @test dataset.tensor_rank == 0
    @test dataset.n_operations == 32
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [1 0 0; 0 0 1; 0 -1 0],
        [-1 0 0; 0 0 -1; 0 1 0],
        [1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 1],
        [1 0 0; 0 0 -1; 0 1 0],
        [-1 0 0; 0 0 1; 0 -1 0],
        [-1 0 0; 0 0 1; 0 1 0],
        [1 0 0; 0 0 -1; 0 -1 0],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
        [-1 0 0; 0 0 -1; 0 -1 0],
        [1 0 0; 0 0 1; 0 1 0],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [1 0 0; 0 0 1; 0 -1 0],
        [-1 0 0; 0 0 -1; 0 1 0],
        [1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 1],
        [1 0 0; 0 0 -1; 0 1 0],
        [-1 0 0; 0 0 1; 0 -1 0],
        [-1 0 0; 0 0 1; 0 1 0],
        [1 0 0; 0 0 -1; 0 -1 0],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
        [-1 0 0; 0 0 -1; 0 -1 0],
        [1 0 0; 0 0 1; 0 1 0],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
    ]  # Compared with Python results
    @test dataset.translations ≈ [
        [0.0, -0.0, 0.0],
        [0.2, 0.2, 0.0],
        [0.0, 0.6, 0.6],
        [0.2, 0.6, 0.4],
        [0.0, 0.2, 0.0],
        [0.2, 0.0, 0.0],
        [0.0, 0.6, 0.4],
        [0.2, 0.6, 0.6],
        [0.2, 0.6, 0.4],
        [0.0, 0.6, 0.6],
        [0.2, 0.2, 0.0],
        [0.0, -0.0, 0.0],
        [0.2, 0.6, 0.6],
        [0.0, 0.6, 0.4],
        [0.2, 0.0, 0.0],
        [0.0, 0.2, 0.0],
        [0.0, 0.5, 0.5],
        [0.2, 0.7, 0.5],
        [0.0, 0.1, 0.1],
        [0.2, 0.1, 0.9],
        [0.0, 0.7, 0.5],
        [0.2, 0.5, 0.5],
        [0.0, 0.1, 0.9],
        [0.2, 0.1, 0.1],
        [0.2, 0.1, 0.9],
        [0.0, 0.1, 0.1],
        [0.2, 0.7, 0.5],
        [0.0, 0.5, 0.5],
        [0.2, 0.1, 0.1],
        [0.0, 0.1, 0.9],
        [0.2, 0.5, 0.5],
        [0.0, 0.7, 0.5],
    ]  # Compared with Python results
    @test dataset.time_reversals == [
        false,
        false,
        true,
        true,
        false,
        false,
        true,
        true,
        true,
        true,
        false,
        false,
        true,
        true,
        false,
        false,
        true,
        true,
        false,
        false,
        true,
        true,
        false,
        false,
        false,
        false,
        true,
        true,
        false,
        false,
        true,
        true,
    ]  # Compared with Python results
    @test dataset.n_atoms == 2
    @test dataset.equivalent_atoms == [0, 0] .+ 1  # Compared with Python results
    @test dataset.transformation_matrix == [
        0.0 0.0 -1.0
        0.0 1.0 0.0
        1.0 0.0 0.0
    ]
    @test dataset.origin_shift == [0.0, 0.9, 0.9]
    @test dataset.n_std_atoms == 2
    @test dataset.std_lattice == Lattice([0.0, 0.0, -4.0], [0.0, 4.0, 0.0], [4.0, 0.0, 0.0])
    @test dataset.std_types == [1, 1]
    @test dataset.std_positions ≈ [[-1.04083409e-17, 0.0, 0.0], [0.5, 0.5, 0.0]]
    @test dataset.std_tensors == [1.0, -1.0]
    @test dataset.std_rotation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.primitive_lattice == Lattice([
        4.0 0.0 0.0
        0.0 4.0 0.0
        0.0 0.0 4.0
    ])
end

# From https://github.com/spglib/spglib/blob/f6abb97/test/functional/fortran/test_fortran_spg_get_symmetry_with_site_tensors.F90#L149-L180
@testset "Test site tensors non-collinear" begin
    lattice = Lattice([
        10 0 0
        0 10 0
        0 0 10
    ])
    positions = [[0.0, 0.0, 0.0]]
    atoms = [1]
    magmoms = [[1, 0, 0]]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 1005
    @test dataset.msg_type == 3
    @test dataset.hall_number == 400
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 16
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 0 1; 0 -1 0],
        [-1 0 0; 0 0 -1; 0 1 0],
        [-1 0 0; 0 0 1; 0 1 0],
        [1 0 0; 0 0 -1; 0 -1 0],
        [-1 0 0; 0 0 -1; 0 -1 0],
        [1 0 0; 0 0 1; 0 1 0],
        [1 0 0; 0 0 -1; 0 1 0],
        [-1 0 0; 0 0 1; 0 -1 0],
    ]
    @test dataset.translations == fill(zeros(3), 16)
    @test dataset.time_reversals == [
        false,
        false,
        true,
        true,
        false,
        false,
        true,
        true,
        false,
        false,
        true,
        true,
        true,
        true,
        false,
        false,
    ]
    @test dataset.n_atoms == 1
    @test dataset.equivalent_atoms == [1]
    @test dataset.transformation_matrix == [
        0.0 0.0 -1.0
        0.0 1.0 0.0
        1.0 0.0 0.0
    ]
    @test dataset.origin_shift == [0.0, 0.0, 0.0]
    @test dataset.n_std_atoms == 1
    @test dataset.std_lattice == Lattice([0, 0, -10], [0, 10, 0], [10, 0, 0])
    @test dataset.std_types == [1]
    @test dataset.std_positions == [[0.0, 0.0, 0.0]]
    @test dataset.std_tensors == [[1.0, 0.0, 0.0]]
    @test dataset.std_rotation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.primitive_lattice == Lattice([
        10 0 0
        0 10 0
        0 0 10
    ])
end
