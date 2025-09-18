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
        @test get_hall_number_from_symmetry(rotations, translations, 1e-5) == 529
    end
    # See https://github.com/spglib/spglib/blob/378240e/python/test/test_collinear_spin.py#L19-L46
    # and https://github.com/spglib/spglib/blob/15a31f6/test/functional/c/test_magnetic_symmetry.cpp#L25-L81
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
            @test equivalent_atoms == [1, 1]
            @testset "Test `get_symmetry_with_site_tensors`" begin
                rotations₂, translations₂, spin_flips = get_symmetry_with_site_tensors(
                    cell, 1e-5
                )
                @test rotations == rotations₂
                @test translations == translations₂
            end
        end
        @testset "Test antiferromagnetism" begin
            magmoms = [1.0, -1.0]
            cell = Cell(lattice, positions, atoms, magmoms)
            rotations, translations, equivalent_atoms = get_symmetry_with_collinear_spin(
                cell, 1e-5
            )
            @test size(rotations) == (96,)
            @test all(iszero(translation) for translation in translations[1:48])
            @test all(
                translation == [1 / 2, 1 / 2, 1 / 2] for translation in translations[49:96]
            )  # Compared with Python
            @test equivalent_atoms == [1, 1]
            @testset "Test `get_symmetry_with_site_tensors`" begin
                rotations₂, translations₂, spin_flips = get_symmetry_with_site_tensors(
                    cell, 1e-5
                )
                @test rotations == rotations₂
                @test translations == translations₂
            end
        end
        @testset "Test broken magmoms" begin
            magmoms = [1.0, 2.0]
            cell = Cell(lattice, positions, atoms, magmoms)
            rotations, translations, equivalent_atoms = get_symmetry_with_collinear_spin(
                cell, 1e-5
            )
            @test size(rotations) == (48,)
            @test size(translations) == (48,)
            @test all(iszero(translation) for translation in translations)  # Compared with Python
            @test equivalent_atoms == [1, 2]  # Compared with Python
            @testset "Test `get_symmetry_with_site_tensors`" begin
                rotations₂, translations₂, spin_flips = get_symmetry_with_site_tensors(
                    cell, 1e-5
                )
                @test rotations == rotations₂
                @test translations == translations₂
            end
        end
    end
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/c/test_magnetic_symmetry.cpp#L84-L160
@testset "Test LaMnO₃" begin
    lattice = Lattice([[5.7461, 0, 0], [0, 7.6637, 0], [0, 0, 5.5333]])
    positions = [
        [0.051300, 0.250000, 0.990500],  # La
        [0.948700, 0.750000, 0.009500],  # La
        [0.551300, 0.250000, 0.509500],  # La
        [0.448700, 0.750000, 0.490500],  # La
        [0.000000, 0.000000, 0.500000],  # Mn
        [0.000000, 0.500000, 0.500000],  # Mn
        [0.500000, 0.500000, 0.000000],  # Mn
        [0.500000, 0.000000, 0.000000],  # Mn
        [0.484900, 0.250000, 0.077700],  # O1
        [0.515100, 0.750000, 0.922300],  # O1
        [0.984900, 0.250000, 0.422300],  # O1
        [0.015100, 0.750000, 0.577700],  # O1
        [0.308500, 0.040800, 0.722700],  # O2
        [0.691500, 0.540800, 0.277300],  # O2
        [0.691500, 0.959200, 0.277300],  # O2
        [0.308500, 0.459200, 0.722700],  # O2
        [0.808500, 0.459200, 0.777300],  # O2
        [0.191500, 0.959200, 0.222700],  # O2
        [0.191500, 0.540800, 0.222700],  # O2
        [0.808500, 0.040800, 0.777300],  # O2
    ]
    atoms = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
    magmoms = [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [3.87, 0, 0],
        [-3.87, 0, 0],
        [-3.87, 0, 0],
        [3.87, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
    ]
    cell = Cell(lattice, positions, atoms, magmoms)
    rotations, translations, spin_flips = get_symmetry_with_site_tensors(cell, 1e-5)
    @test length(rotations) == length(translations) == length(spin_flips) == 8
    @test rotations == [
        [
            1 0 0
            0 1 0
            0 0 1
        ],
        [
            -1 0 0
            0 -1 0
            0 0 -1
        ],
        [
            -1 0 0
            0 -1 0
            0 0 1
        ],
        [
            1 0 0
            0 1 0
            0 0 -1
        ],
        [
            1 0 0
            0 -1 0
            0 0 -1
        ],
        [
            -1 0 0
            0 1 0
            0 0 1
        ],
        [
            -1 0 0
            0 1 0
            0 0 -1
        ],
        [
            1 0 0
            0 -1 0
            0 0 1
        ],
    ]  # Compared with Python results
    @test translations == [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.0, 0.5, 0.0],
        [0.0, 0.5, 0.0],
    ]  # Compared with Python results
    dataset = get_magnetic_dataset(cell, 1e-5)
    @test rotations == dataset.rotations
    @test translations == dataset.translations
    @test dataset.time_reversals == [false, false, true, true, true, true, false, false]  # Compared with Python results
    @test dataset.time_reversals == (1 .- spin_flips) ./ 2  # See https://github.com/spglib/spglib/blob/v2.1.0/test/functional/c/test_magnetic_symmetry.cpp#L156-L160
    @testset "Test `get_magnetic_spacegroup_type_from_symmetry`" begin
        # From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/c/test_magnetic_symmetry.cpp#L161-L164
        spgtype = get_magnetic_spacegroup_type_from_symmetry(
            rotations, translations, dataset.time_reversals, lattice, 1e-5
        )
        @test spgtype.uni_number == 546
    end
end

# See https://github.com/singularitti/Spglib.jl/issues/91#issuecomment-1206106977
@testset "Test example given by Jae-Mo Lihm (@jaemolihm)" begin
    lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    positions = [[-0.1, -0.1, -0.1], [0.1, 0.1, 0.1]]
    atoms = [1, 1]
    magmoms = [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]
    cell = Cell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell, 1e-5)
    @test dataset.uni_number == 68
    @test dataset.msg_type == 3
    @test dataset.hall_number == 63
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 4
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [0 -1 0; -1 0 0; 0 0 -1],
        [0 1 0; 1 0 0; 0 0 1],
    ]  # Compared with Python results
    @test dataset.translations == fill(zeros(3), 4)  # Compared with Python results
    @test dataset.time_reversals == [false, true, false, true]  # Compared with Python results
    @test dataset.n_atoms == 2
    @test dataset.equivalent_atoms == [0, 0] .+ 1  # Compared with Python results
    @test dataset.transformation_matrix == [
        0.5 0.5 0.0
        -0.5 0.5 0.0
        0.0 0.0 1.0
    ]
    @test dataset.origin_shift == [0.0, 0.0, 0.0]
    @test dataset.n_std_atoms == 4
    @test dataset.std_lattice == Lattice([
        1.0 -1.0 0.0
        1.0 1.0 0.0
        0.0 0.0 1.0
    ])
    @test dataset.std_types == [1, 1, 1, 1]
    @test dataset.std_positions ≈
        [[0.9, 0.0, 0.9], [0.4, 0.5, 0.9], [0.1, 0.0, 0.1], [0.6, 0.5, 0.1]]
    @test dataset.std_tensors ==
        [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [0.0, 0.0, -1.0]]
    @test dataset.std_rotation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.primitive_lattice == Lattice([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_magnetic_dataset.py#L9-L44
@testset "Test type I" begin
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
    @test dataset.origin_shift ≈ [5.55111512e-17, -1.66533454e-16, 0]
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
    @test dataset.std_tensors ≈ [
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

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_magnetic_dataset.py#L46-L87
@testset "Test type III" begin
    lattice = [
        [1.06949000e01, 0.00000000e00, 0.00000000e00],
        [3.84998337e-16, 6.28750000e00, 0.00000000e00],
        [3.09590711e-16, 3.09590711e-16, 5.05600000e00],
    ]
    positions = [
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.5],
        [0.0, 0.5, 0.0],
        [0.2794, 0.25, 0.9906],
        [0.2206, 0.75, 0.4906],
        [0.7206, 0.75, 0.0094],
        [0.7794, 0.25, 0.5094],
    ]
    atoms = [0, 0, 0, 0, 0, 0, 0, 0]
    magmoms = [
        [3.0, 0.4, 0.0],
        [-3.0, -0.4, 0.0],
        [-3.0, 0.4, 0.0],
        [3.0, -0.4, 0.0],
        [4.5, 0.0, 0.0],
        [-4.5, 0.0, 0.0],
        [4.5, 0.0, 0.0],
        [-4.5, 0.0, 0.0],
    ]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 544
    @test dataset.msg_type == 3
    @test dataset.hall_number == 292
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 8
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
    ]  # Compared with Python results
    @test dataset.translations ≈ [
        [0.0, 0.0, 0.0],
        [-1.06253018e-65, -9.77423530e-34, 0.0],
        [0.5, -1.46613529e-33, 0.5],
        [0.5, -4.88711765e-34, 0.5],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [5.43536217e-33, 0.5, 0.0],
        [5.43536217e-33, 0.5, 0.0],
    ]  # Compared with Python results
    @test dataset.time_reversals == [false, false, false, false, true, true, true, true]  # Compared with Python results
    @test dataset.n_atoms == 8
    @test dataset.equivalent_atoms == [0, 0, 0, 0, 4, 4, 4, 4] .+ 1  # Compared with Python results
    @test dataset.transformation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.origin_shift ≈ [0.0, 7.33067647e-34, 0.0] atol = eps()
    @test dataset.n_std_atoms == 8
    @test dataset.std_lattice ≈ Lattice([
        [10.6949000, 0.0, 0.0],
        [3.84998337e-16, 6.2875, 0.0],
        [3.09590711e-16, 3.09590711e-16, 5.056],
    ])
    @test dataset.std_types == [1, 1, 1, 1, 1, 1, 1, 1]  # Python results are all zeros?
    @test dataset.std_positions ≈ [
        [1.35884054e-33, 3.66533824e-34, 0.0],
        [0.5, 1.34395735e-33, 0.5],
        [0.5, 0.5, 0.5],
        [-4.90185411e-33, 0.5, 0.0],
        [0.2794, 0.25, 0.9906],
        [0.2206, 0.75, 0.4906],
        [0.7206, 0.75, 0.0094],
        [0.7794, 0.25, 0.5094],
    ]
    @test dataset.std_tensors == [
        [3.0, 0.4, 0.0],
        [-3.0, -0.4, 0.0],
        [-3.0, 0.4, 0.0],
        [3.0, -0.4, 0.0],
        [4.5, 0.0, 0.0],
        [-4.5, 0.0, 0.0],
        [4.5, 0.0, 0.0],
        [-4.5, 0.0, 0.0],
    ]
    @test dataset.std_rotation_matrix ≈ [
        1.0 9.31289551e-33 1.64960847e-35
        0.0 1.0 1.64960847e-35
        0.0 0.0 1.0
    ]
    @test dataset.primitive_lattice == Lattice([
        [-3.09590711e-16, -3.09590711e-16, -5.056],
        [-3.84998337e-16, -6.2875, 0],
        [-1.06949000e+01, 0, 0],
    ])
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_magnetic_dataset.py#L89-L146
@testset "Test monoclinic" begin
    a = 10.6926
    b = 6.2851
    c = 5.0557
    lattice = Lattice([10.6926, 0.0, 0.0], [0.0, 6.2851, 0.0], [0.0, 0.0, 5.0557])
    positions = [
        [0, 0, 0],
        [0.5, 0, 0.5],
        [0, 0.5, 0],
        [0.5, 0.5, 0.5],
        [0.2794, 0.25, 0.9903],
        [0.2206, 0.75, 0.4903],
        [0.7206, 0.75, 0.0097],
        [0.7794, 0.25, 0.5097],
    ]
    atoms = [0, 0, 0, 0, 0, 0, 0, 0]
    magmoms = [
        [1.9, 2.8, 0.3],
        [-1.9, -2.8, 0.3],
        [1.7, 2.4, -0.3],
        [-1.7, -2.4, -0.3],
        [2.9, 3.4, 0],
        [-2.9, -3.4, 0],
        [2.9, 3.4, 0],
        [-2.9, -3.4, 0],
    ]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 82
    @test dataset.msg_type == 1
    @test dataset.hall_number == 81
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 4
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
    ]  # Compared with Python results
    @test dataset.translations ==
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.5, 0.0, 0.5], [0.5, 0.0, 0.5]]  # Compared with Python results
    @test dataset.time_reversals == [false, false, false, false]  # Compared with Python results
    @test dataset.n_atoms == 8
    @test dataset.equivalent_atoms == [0, 0, 2, 2, 4, 4, 4, 4] .+ 1  # Compared with Python results
    @test dataset.transformation_matrix == [
        0.0 1.0 0.0
        0.0 0.0 1.0
        1.0 0.0 0.0
    ]
    @test dataset.origin_shift == [0.0, 0.0, 0.0]
    @test dataset.n_std_atoms == 8
    @test dataset.std_lattice ==
        Lattice([[0.0, 6.2851, 0.0], [0.0, 0.0, 5.0557], [10.6926, 0.0, 0.0]])
    @test dataset.std_types == [1, 1, 1, 1, 1, 1, 1, 1]  # Python results are all zeros?
    @test dataset.std_positions ≈ [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.25, 0.9903, 0.2794],
        [0.75, 0.4903, 0.2206],
        [0.75, 0.0097, 0.7206],
        [0.25, 0.5097, 0.7794],
    ]
    @test dataset.std_tensors == [
        [1.9, 2.8, 0.3],
        [-1.9, -2.8, 0.3],
        [1.7, 2.4, -0.3],
        [-1.7, -2.4, -0.3],
        [2.9, 3.4, 0.0],
        [-2.9, -3.4, 0.0],
        [2.9, 3.4, 0.0],
        [-2.9, -3.4, 0.0],
    ]
    @test dataset.std_rotation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.primitive_lattice ==
        Lattice([[0.0, 0.0, -5.0557], [0.0, -6.2851, 0.0], [-10.6926, 0.0, 0.0]])
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_magnetic_dataset.py#L148-L200
@testset "Test centered monoclinic" begin
    lattice = [[8.627, 0.0, 0.0], [0.0, 8.627, 0.0], [0.0, 0.0, 13.822]]
    positions = [
        [0.0925, 0.19644, 0.24573],
        [0.5925, 0.69644, 0.74573],
        [0.0925, 0.19644, 0.75427],
        [0.5925, 0.69644, 0.25427],
        [0.9075, 0.80356, 0.24573],
        [0.4075, 0.30356, 0.74573],
        [0.9075, 0.80356, 0.75427],
        [0.4075, 0.30356, 0.25427],
        [0.80356, 0.0925, 0.24573],
        [0.30356, 0.5925, 0.74573],
        [0.80356, 0.0925, 0.75427],
        [0.30356, 0.5925, 0.25427],
        [0.19644, 0.9075, 0.24573],
        [0.69644, 0.4075, 0.74573],
        [0.19644, 0.9075, 0.75427],
        [0.69644, 0.4075, 0.25427],
    ]
    atoms = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    magmoms = [
        [2.674, -1.337, 0.0],
        [2.674, -1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [2.674, -1.337, 0.0],
        [2.674, -1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [2.674, -1.337, 0.0],
        [2.674, -1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [2.674, -1.337, 0.0],
        [2.674, -1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [-2.674, 1.337, 0.0],
    ]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 67
    @test dataset.msg_type == 3
    @test dataset.hall_number == 63
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 8
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
    ]
    @test dataset.translations ≈ [
        [0.0, 0.0, 0.0],
        [5.55111512e-17, -5.55111512e-17, -5.55111512e-17],
        [5.55111512e-17, -5.55111512e-17, 0.0],
        [0.0, 0.0, -5.55111512e-17],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
    ]
    @test dataset.time_reversals == [false, true, true, false, false, true, true, false]
    @test dataset.n_atoms == 16
    @test dataset.equivalent_atoms == [0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8] .+ 1
    @test dataset.transformation_matrix == [
        1.0 0.0 0.0
        0.0 0.0 -1.0
        1.0 1.0 0.0
    ]
    @test dataset.origin_shift ≈ [-2.77555756e-17, 2.77555756e-17, 6.16297582e-33]
    @test dataset.n_std_atoms == 16
    @test dataset.std_lattice ≈
        Lattice([[8.627, -8.627, 0.0], [0.0, 0.0, -13.822], [0.0, 8.627, 0.0]])
    @test dataset.std_types == fill(1, 16)  # Python use all zeros?
    @test dataset.std_positions ≈ [
        [0.0925, 0.75427, 0.28894],
        [0.5925, 0.25427, 0.28894],
        [0.5925, 0.74573, 0.28894],
        [0.0925, 0.24573, 0.28894],
        [0.9075, 0.75427, 0.71106],
        [0.4075, 0.25427, 0.71106],
        [0.4075, 0.74573, 0.71106],
        [0.9075, 0.24573, 0.71106],
        [0.80356, 0.75427, 0.89606],
        [0.30356, 0.25427, 0.89606],
        [0.30356, 0.74573, 0.89606],
        [0.80356, 0.24573, 0.89606],
        [0.19644, 0.75427, 0.10394],
        [0.69644, 0.25427, 0.10394],
        [0.69644, 0.74573, 0.10394],
        [0.19644, 0.24573, 0.10394],
    ]
    @test dataset.std_tensors ≈ [
        [2.674, -1.337, 0.0],
        [2.674, -1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [2.674, -1.337, 0.0],
        [2.674, -1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [2.674, -1.337, 0.0],
        [2.674, -1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [2.674, -1.337, 0.0],
        [2.674, -1.337, 0.0],
        [-2.674, 1.337, 0.0],
        [-2.674, 1.337, 0.0],
    ]
    @test dataset.std_rotation_matrix ≈ [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.primitive_lattice ≈
        Lattice([[-8.627, -0.0, -0.0], [-0.0, -8.627, -0.0], [4.3135, 4.3135, 6.911]])
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_magnetic_dataset.py#L202C9-L240
@testset "Test trigonal" begin
    lattice = [[3.77, 0.0, 0.0], [0.0, 3.77, 0.0], [0.0, 0.0, 3.77]]
    positions = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
    atoms = [0, 1, 1, 1]
    magmoms = [[0.0, 0.0, 0.0], [2.0, -1.0, -1.0], [-1.0, 2.0, -1.0], [-1.0, -1.0, 2.0]]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 1331
    @test dataset.msg_type == 3
    @test dataset.hall_number == 458
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 12
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [0 -1 0; -1 0 0; 0 0 -1],
        [0 1 0; 1 0 0; 0 0 1],
        [0 0 1; 1 0 0; 0 1 0],
        [0 0 -1; -1 0 0; 0 -1 0],
        [0 0 -1; 0 -1 0; -1 0 0],
        [0 0 1; 0 1 0; 1 0 0],
        [0 1 0; 0 0 1; 1 0 0],
        [0 -1 0; 0 0 -1; -1 0 0],
        [-1 0 0; 0 0 -1; 0 -1 0],
        [1 0 0; 0 0 1; 0 1 0],
    ]  # Compared with Python results
    @test dataset.translations == [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]  # Compared with Python results
    @test dataset.time_reversals ==
        [false, false, true, true, false, false, true, true, false, false, true, true]  # Compared with Python results
    @test dataset.n_atoms == 4
    @test dataset.equivalent_atoms == [0, 1, 1, 1] .+ 1  # Compared with Python results
    @test dataset.transformation_matrix ≈ [
        2 -1 -1
        1 1 -2
        1 1 1
    ] / 3
    @test dataset.origin_shift == [0.0, 0.0, 0.0]
    @test dataset.n_std_atoms == 12
    @test dataset.std_lattice ≈ Lattice([
        [3.77, -3.77, -6.97590134e-17], [6.97590134e-17, 3.77, -3.77], [3.77, 3.77, 3.77]
    ])
    @test dataset.std_types == [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] .+ 1  # Python results are different
    @test dataset.std_positions ≈ [
        [0.0, 0.0, 0.0],
        [0.66666667, 0.33333333, 0.33333333],
        [0.33333333, 0.66666667, 0.66666667],
        [0.66666667, 0.83333333, 0.33333333],
        [0.33333333, 0.16666667, 0.66666667],
        [5.55111512e-17, 0.5, 0.0],
        [0.16666667, 0.83333333, 0.33333333],
        [0.83333333, 0.16666667, 0.66666667],
        [0.5, 0.5, 0.0],
        [0.16666667, 0.33333333, 0.33333333],
        [0.83333333, 0.66666667, 0.66666667],
        [0.5, 0.0, 0.0],
    ]
    @test dataset.std_tensors ≈ [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [2.0, -1.0, -1.0],
        [2.0, -1.0, -1.0],
        [2.0, -1.0, -1.0],
        [-1.0, 2.0, -1.0],
        [-1.0, 2.0, -1.0],
        [-1.0, 2.0, -1.0],
        [-1.0, -1.0, 2.0],
        [-1.0, -1.0, 2.0],
        [-1.0, -1.0, 2.0],
    ]
    @test dataset.std_rotation_matrix ≈ [
        1.0 1.12887399e-18 1.12887399e-18
        1.12887399e-18 1.0 1.12887399e-18
        1.12887399e-18 1.12887399e-18 1.0
    ]
    @test dataset.primitive_lattice == Lattice([
        3.77 0.0 0.0
        0.0 3.77 0.0
        0.0 0.0 3.77
    ])
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_magnetic_dataset.py#L257-L319
@testset "Test conventional" begin
    lattice = [[4.957, 0.0, 0.0], [-2.4785, 4.29288793, 0.0], [0.0, 0.0, 13.5923]]
    positions = [
        [0.0, 0.0, 0.34751],
        [1 / 3, 2 / 3, 0.01417667],
        [2 / 3, 1 / 3, 0.68084333],
        [1 / 3, 2 / 3, 0.51417667],
        [2 / 3, 1 / 3, 0.18084333],
        [0.0, 0.0, 0.84751],
        [1 / 3, 2 / 3, 0.81915667],
        [0.0, 0.0, 0.15249],
        [2 / 3, 1 / 3, 0.48582333],
        [0.0, 0.0, 0.65249],
        [2 / 3, 1 / 3, 0.98582333],
        [1 / 3, 2 / 3, 0.31915667],
        [0.3057, 0.0, 0.25],
        [0.63903333, 2 / 3, 0.91666667],
        [0.97236667, 1 / 3, 0.58333333],
        [0.02763333, 2 / 3, 0.41666667],
        [0.36096667, 1 / 3, 0.08333333],
        [0.6943, 0.0, 0.75],
    ]
    atoms = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
    magmoms = [
        [-1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 94
    @test dataset.msg_type == 3
    @test dataset.hall_number == 90
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 12
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [1 -1 0; 0 -1 0; 0 0 -1],
        [-1 1 0; 0 1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [1 -1 0; 0 -1 0; 0 0 -1],
        [-1 1 0; 0 1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [1 -1 0; 0 -1 0; 0 0 -1],
        [-1 1 0; 0 1 0; 0 0 1],
    ]  # Compared with Python results
    @test dataset.translations ≈ [
        [0.0, 0.0, 0.0],
        [0.0, -1.11022302e-16, 0.0],
        [0.66666667, 0.33333333, 0.83333333],
        [0.66666667, 0.33333333, 0.83333333],
        [0.33333333, 0.66666667, 0.66666667],
        [0.33333333, 0.66666667, 0.66666667],
        [0.0, 0.0, 0.5],
        [-5.55111512e-17, 0.0, 0.5],
        [0.66666667, 0.33333333, 0.33333333],
        [0.66666667, 0.33333333, 0.33333333],
        [0.33333333, 0.66666667, 0.16666667],
        [0.33333333, 0.66666667, 0.16666667],
    ]  # Compared with Python results
    @test dataset.time_reversals ==
        [false, true, true, false, false, true, true, false, false, true, true, false]  # Compared with Python results
    @test dataset.n_atoms == 18
    @test dataset.equivalent_atoms ==
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 12, 12, 12, 12, 12] .+ 1
    @test dataset.transformation_matrix ≈ [
        8.32667268e-17 0.5 1.0
        -1.0 0.5 -9.25185854e-17
        0.0 -1.0 1.0
    ]
    @test dataset.origin_shift ≈ [0.25, 0.75, -5.55111512e-17]
    @test dataset.n_std_atoms == 12
    @test dataset.std_lattice ≈ Lattice([
        [-2.50959747e-16, 2.86192529, 9.06153333],
        [-4.957, 2.38303151e-16, 7.54524221e-16],
        [-1.78604045e-16, -2.86192529, 4.53076667],
    ])
    @test dataset.std_types == [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1] .+ 1
    @test dataset.std_positions ≈ [
        [0.59751, 0.75, 0.34751],
        [0.09751, 0.25, 0.34751],
        [0.59751, 0.25, 0.84751],
        [0.09751, 0.75, 0.84751],
        [0.90249, 0.25, 0.15249],
        [0.40249, 0.75, 0.15249],
        [0.90249, 0.75, 0.65249],
        [0.40249, 0.25, 0.65249],
        [0.0, 0.9443, 0.25],
        [0.5, 0.4443, 0.25],
        [0.5, 0.5557, 0.75],
        [0.0, 0.0557, 0.75],
    ]
    @test dataset.std_tensors ≈ [
        [-1.0, -3.43333113e-34, -7.30671738e-33],
        [-1.0, -3.43333113e-34, -7.30671738e-33],
        [-1.0, -3.43333113e-34, -7.30671738e-33],
        [-1.0, -3.43333113e-34, -7.30671738e-33],
        [1.0, 3.43333113e-34, 7.30671738e-33],
        [1.0, 3.43333113e-34, 7.30671738e-33],
        [1.0, 3.43333113e-34, 7.30671738e-33],
        [1.0, 3.43333113e-34, 7.30671738e-33],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
    @test dataset.std_rotation_matrix ≈ [
        1.0 -4.88817852e-17 1.06865048e-18
        3.43333113e-34 1.0 0.0
        7.30671738e-33 2.06895319e-16 1.0
    ]
    @test dataset.primitive_lattice ≈ Lattice([
        [4.957, -0.0, -0.0],
        [-2.4785, -4.29288793, -0.0],
        [2.4785, -1.43096264, -4.53076667],
    ])
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_magnetic_dataset.py#L321-L368
@testset "Test type IV" begin
    lattice = [[16.0831, 0.0, 0.0], [0.0, 4.3887, 0.0], [-5.97205829, 0.0, 12.38135828]]
    positions = [
        [0.0, 0.0, 0.25],
        [0.5, 0.5, 0.25],
        [0.0, 0.0, 0.75],
        [0.5, 0.5, 0.75],
        [0.129, 0.0, 0.0585],
        [0.629, 0.5, 0.0585],
        [0.129, 0.0, 0.5585],
        [0.629, 0.5, 0.5585],
        [0.871, 0.0, 0.4415],
        [0.371, 0.5, 0.4415],
        [0.871, 0.0, 0.9415],
        [0.371, 0.5, 0.9415],
    ]
    atoms = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    magmoms = [
        [0.0, 8.86, 0.0],
        [0.0, 8.86, 0.0],
        [0.0, -8.86, 0.0],
        [0.0, -8.86, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 70
    @test dataset.msg_type == 4
    @test dataset.hall_number == 63
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 16
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
    ]  # Compared with Python results
    @test dataset.translations == [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.0, 0.0, 0.5],
        [0.0, 0.0, 0.5],
        [0.0, 0.0, 0.5],
        [0.0, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.0],
    ]  # Compared with Python results
    @test dataset.time_reversals == [
        false,
        true,
        true,
        false,
        true,
        false,
        false,
        true,
        true,
        false,
        false,
        true,
        false,
        true,
        true,
        false,
    ]  # Compared with Python results
    @test dataset.n_atoms == 12
    @test dataset.equivalent_atoms == [0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4] .+ 1  # Compared with Python results
    @test dataset.transformation_matrix == [
        -1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 -1.0
    ]
    @test dataset.origin_shift == [0.0, 0.0, 0.75]
    @test dataset.n_std_atoms == 12
    @test dataset.std_lattice ≈ Lattice([
        [-16.0831, 0.0, 0.0], [0.0, 4.3887, 0.0], [5.97205829, 0.0, -12.38135828]
    ])
    @test dataset.std_types == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    @test dataset.std_positions ≈ [
        [1.38777878e-17, 1.38777878e-17, 0.5],
        [0.5, 0.5, 0.5],
        [0.0, 0.0, 1.11022302e-16],
        [0.5, 0.5, 2.22044605e-16],
        [0.371, 0.5, 0.6915],
        [0.871, 0.0, 0.6915],
        [0.871, 0.0, 0.1915],
        [0.371, 0.5, 0.1915],
        [0.129, 2.77555756e-17, 0.3085],
        [0.629, 0.5, 0.3085],
        [0.629, 0.5, 0.8085],
        [0.129, 0.0, 0.8085],
    ]
    @test dataset.std_tensors ≈ [
        [0.0, 8.86, 0.0],
        [0.0, 8.86, 0.0],
        [0.0, -8.86, 0.0],
        [0.0, -8.86, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
    @test dataset.std_rotation_matrix ≈ [
        1.0 0.0 1.32046516e-16
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]
    @test dataset.primitive_lattice ≈ Lattice([
        [0.0, 4.3887, 0.0], [-8.04155, -2.19435, 0.0], [2.06949171, -2.19435, 12.38135828]
    ])
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_magnetic_dataset.py#L370-L412
@testset "Test nonstandard setting" begin
    lattice = [[6.1096, 0.0, 0.0], [0.0, 7.2555, 0.0], [0.0, 0.0, 7.5708]]
    positions = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.2732, 0.8579, 0.352],
        [0.2268, 0.3579, 0.148],
        [0.7268, 0.1421, 0.648],
        [0.7732, 0.6421, 0.852],
        [0.7732, 0.6421, 0.148],
        [0.7268, 0.1421, 0.352],
        [0.2268, 0.3579, 0.852],
        [0.2732, 0.8579, 0.648],
    ]
    atoms = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
    magmoms = [
        [5.01, 0.0, 0.0],
        [-5.01, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 496
    @test dataset.msg_type == 3
    @test dataset.hall_number == 275
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 8
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
    ]  # Compared with Python results
    @test dataset.translations ≈ [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, -1.11022302e-16],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, -1.11022302e-16],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
    ]  # Compared with Python results
    @test dataset.time_reversals == [false, false, true, true, true, true, false, false]  # Compared with Python results
    @test dataset.n_atoms == 10
    @test dataset.equivalent_atoms == [0, 0, 2, 2, 2, 2, 2, 2, 2, 2] .+ 1  # Compared with Python results
    @test dataset.transformation_matrix == [
        0.0 -1.0 0.0
        -1.0 0.0 0.0
        0.0 0.0 -1.0
    ]
    @test dataset.origin_shift == [0.0, 0.0, 0.0]
    @test dataset.n_std_atoms == 10
    @test dataset.std_lattice ≈
        Lattice([[-0.0, -7.2555, -0.0], [-6.1096, -0.0, -0.0], [-0.0, -0.0, -7.5708]])
    @test dataset.std_types == [0, 0, 1, 1, 1, 1, 1, 1, 1, 1] .+ 1
    @test dataset.std_positions ≈ [
        [-2.77555756e-17, 0.0, 6.16297582e-33],
        [0.5, 0.5, 0.5],
        [0.1421, 0.7268, 0.648],
        [0.6421, 0.7732, 0.852],
        [0.8579, 0.2732, 0.352],
        [0.3579, 0.2268, 0.148],
        [0.3579, 0.2268, 0.852],
        [0.8579, 0.2732, 0.648],
        [0.6421, 0.7732, 0.148],
        [0.1421, 0.7268, 0.352],
    ]
    @test dataset.std_tensors ≈ [
        [5.01, 0.0, 0.0],
        [-5.01, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
    @test dataset.std_rotation_matrix ≈ [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.primitive_lattice ==
        Lattice([[6.1096, 0.0, 0.0], [0.0, 7.2555, 0.0], [0.0, 0.0, 7.5708]])
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_magnetic_dataset.py#L414-L456
@testset "Test nonstandard setting 2" begin
    lattice = [[7.6476, 0.0, 0.0], [0.0, 6.8855, 0.0], [0.0, 0.0, 14.6058]]
    positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.0],
        [0.5, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.3386, 0.25, 0.6426],
        [0.6614, 0.75, 0.3574],
        [0.8386, 0.25, 0.8574],
        [0.1614, 0.75, 0.1426],
        [0.2055, 0.25, 0.4372],
        [0.7945, 0.75, 0.5628],
        [0.7055, 0.25, 0.0628],
        [0.2945, 0.75, 0.9372],
    ]
    atoms = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    magmoms = [
        [0.7, 0.0, 0.7],
        [-0.7, -0.0, -0.7],
        [-0.7, 0.0, 0.7],
        [0.7, -0.0, -0.7],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.7, 0.0, 0.7],
        [-0.7, -0.0, -0.7],
        [-0.7, 0.0, 0.7],
        [0.7, -0.0, -0.7],
    ]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 125
    @test dataset.msg_type == 3
    @test dataset.hall_number == 115
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 4
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
    ]  # Compared with Python results
    @test dataset.translations ≈
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.5], [0.5, 0.5, 0.5], [0.0, 0.5, 0.0]]  # Compared with Python results
    @test dataset.time_reversals == [false, true, true, false]  # Compared with Python results
    @test dataset.n_atoms == 12
    @test dataset.equivalent_atoms == [0, 0, 0, 0, 4, 4, 4, 4, 8, 8, 8, 8] .+ 1  # Compared with Python results
    @test dataset.transformation_matrix == [
        -1.0 0.0 0.0
        0.0 0.0 -1.0
        0.0 -1.0 0.0
    ]
    @test dataset.origin_shift == [0.25, -0.5, 0.25]
    @test dataset.n_std_atoms == 12
    @test dataset.std_lattice ≈
        Lattice([[-7.6476, 0.0, 0.0], [0.0, 0.0, -14.6058], [0.0, -6.8855, 0.0]])
    @test dataset.std_types == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] .+ 1
    @test dataset.std_positions ≈ [
        [0.25, 0.5, 0.25],
        [0.25, 0.5, 0.75],
        [0.75, -5.55111512e-17, 0.75],
        [0.75, -5.55111512e-17, 0.25],
        [0.9114, 0.8574, 0.0],
        [0.5886, 0.1426, 0.5],
        [0.4114, 0.6426, 0.0],
        [0.0886, 0.3574, 0.5],
        [0.0445, 0.0628, 0.0],
        [0.4555, 0.9372, 0.5],
        [0.5445, 0.4372, 0.0],
        [0.9555, 0.5628, 0.5],
    ]
    @test dataset.std_tensors ≈ [
        [0.7, 0.0, 0.7],
        [-0.7, 0.0, -0.7],
        [-0.7, 0.0, 0.7],
        [0.7, 0.0, -0.7],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.7, 0.0, 0.7],
        [-0.7, 0.0, -0.7],
        [-0.7, 0.0, 0.7],
        [0.7, 0.0, -0.7],
    ]
    @test dataset.std_rotation_matrix ≈ [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.primitive_lattice ==
        Lattice([[0, -6.8855, 0], [-7.6476, 0, 0], [0, 0, -14.6058]])
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_magnetic_dataset.py#L461-L512
@testset "Test nonstandard setting 3" begin
    lattice = [[19.2776, 0.0, 0.0], [0.0, 9.2409, 0.0], [0.0, 0.0, 5.4676]]
    positions = [
        [0.1228, 0.3456, 0.3592],
        [0.3772, 0.6544, 0.8592],
        [0.3772, 0.8456, 0.3592],
        [0.1228, 0.1544, 0.8592],
        [0.6228, 0.1544, 0.6408],
        [0.8772, 0.8456, 0.1408],
        [0.8772, 0.6544, 0.6408],
        [0.6228, 0.3456, 0.1408],
        [0.3777, 0.4818, 0.3487],
        [0.1223, 0.5182, 0.8487],
        [0.1223, 0.9818, 0.3487],
        [0.3777, 0.0182, 0.8487],
        [0.8777, 0.0182, 0.6513],
        [0.6223, 0.9818, 0.1513],
        [0.6223, 0.5182, 0.6513],
        [0.8777, 0.4818, 0.1513],
    ]
    atoms = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    magmoms = [
        [4.05, 0, 0.1],
        [-4.05, 0, 0.1],
        [4.05, 0, -0.1],
        [-4.05, 0, -0.1],
        [-4.05, 0, 0.1],
        [4.05, 0, 0.1],
        [-4.05, 0, -0.1],
        [4.05, 0, -0.1],
        [4.05, 0, 0.1],
        [-4.05, 0, 0.1],
        [4.05, 0, -0.1],
        [-4.05, 0, -0.1],
        [-4.05, 0, 0.1],
        [4.05, 0, 0.1],
        [-4.05, 0, -0.1],
        [4.05, 0, -0.1],
    ]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 533
    @test dataset.msg_type == 3
    @test dataset.hall_number == 290
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 8
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 1],
    ]  # Compared with Python results
    @test dataset.translations == [
        [-0.0, -0.0, -0.0],
        [-0.0, -0.0, -0.0],
        [-0.0, 0.5, 0.5],
        [-0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.0],
    ]  # Compared with Python results
    @test dataset.time_reversals == [false, true, true, false, false, true, true, false]  # Compared with Python results
    @test dataset.n_atoms == 16
    @test dataset.equivalent_atoms == [0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8] .+ 1  # Compared with Python results
    @test dataset.transformation_matrix == [0.0 0.0 -1.0; -1.0 0.0 0.0; 0.0 1.0 0.0]
    @test dataset.origin_shift == [0.0, 0.0, 0.0]
    @test dataset.n_std_atoms == 16
    @test dataset.std_lattice ≈
        Lattice([[0.0, 0.0, -5.4676], [-19.2776, 0.0, 0.0], [0.0, 9.2409, 0.0]])  # To work on 32-bit systems
    @test dataset.std_types == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    @test dataset.std_positions ≈ [
        [0.6408, 0.8772, 0.3456],
        [0.1408, 0.6228, 0.6544],
        [0.6408, 0.6228, 0.8456],
        [0.1408, 0.8772, 0.1544],
        [0.3592, 0.3772, 0.1544],
        [0.8592, 0.1228, 0.8456],
        [0.3592, 0.1228, 0.6544],
        [0.8592, 0.3772, 0.3456],
        [0.6513, 0.6223, 0.4818],
        [0.1513, 0.8777, 0.5182],
        [0.6513, 0.8777, 0.9818],
        [0.1513, 0.6223, 0.0182],
        [0.3487, 0.1223, 0.0182],
        [0.8487, 0.3777, 0.9818],
        [0.3487, 0.3777, 0.5182],
        [0.8487, 0.1223, 0.4818],
    ]
    @test dataset.std_tensors ≈ [
        [4.05, 0.0, 0.1],
        [-4.05, 0.0, 0.1],
        [4.05, 0.0, -0.1],
        [-4.05, 0.0, -0.1],
        [-4.05, 0.0, 0.1],
        [4.05, 0.0, 0.1],
        [-4.05, 0.0, -0.1],
        [4.05, 0.0, -0.1],
        [4.05, 0.0, 0.1],
        [-4.05, 0.0, 0.1],
        [4.05, 0.0, -0.1],
        [-4.05, 0.0, -0.1],
        [-4.05, 0.0, 0.1],
        [4.05, 0.0, 0.1],
        [-4.05, 0.0, -0.1],
        [4.05, 0.0, -0.1],
    ]
    @test dataset.std_rotation_matrix ≈ [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test dataset.primitive_lattice ==
        Lattice([[0.0, 0.0, -5.4676], [0.0, -9.2409, 0.0], [-19.2776, 0.0, 0.0]])
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_magnetic_dataset.py#L514-L585
@testset "Test nonstandard type 4" begin
    lattice = [[7.5091, 0.0, 0.0], [0.0, 5.7319, 0.0], [0.0, 0.0, 22.543]]
    positions = [
        [0.25, 0.0, 0.10125],
        [0.5, 0.5, 0.85125],
        [0.75, 0.0, 0.60125],
        [0.0, 0.5, 0.35125],
        [0.75, 0.0, 0.10125],
        [0.0, 0.5, 0.85125],
        [0.5, 0.5, 0.35125],
        [0.25, 0.0, 0.60125],
        [0.75, 0.0, 0.39875],
        [0.5, 0.5, 0.64875],
        [0.25, 0.0, 0.89875],
        [0.0, 0.5, 0.14875],
        [0.25, 0.0, 0.39875],
        [0.0, 0.5, 0.64875],
        [0.5, 0.5, 0.14875],
        [0.75, 0.0, 0.89875],
        [0.0, 0.0, 0.0],
        [0.25, 0.5, 0.75],
        [0.5, 0.0, 0.5],
        [0.75, 0.5, 0.25],
        [0.5, 0.0, 0.0],
        [0.75, 0.5, 0.75],
        [0.25, 0.5, 0.25],
        [0.0, 0.0, 0.5],
    ]
    atoms = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
    magmoms = [
        [0.12, 0, 9.06],
        [0.12, 0, 9.06],
        [0.12, 0, 9.06],
        [0.12, 0, 9.06],
        [-0.12, 0, -9.06],
        [-0.12, 0, -9.06],
        [-0.12, 0, -9.06],
        [-0.12, 0, -9.06],
        [-0.12, 0, -9.06],
        [-0.12, 0, -9.06],
        [-0.12, 0, -9.06],
        [-0.12, 0, -9.06],
        [0.12, 0, 9.06],
        [0.12, 0, 9.06],
        [0.12, 0, 9.06],
        [0.12, 0, 9.06],
        [0.58, 0, -1.26],
        [0.58, 0, -1.26],
        [0.58, 0, -1.26],
        [0.58, 0, -1.26],
        [-0.58, 0, 1.26],
        [-0.58, 0, 1.26],
        [-0.58, 0, 1.26],
        [-0.58, 0, 1.26],
    ]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    dataset = get_magnetic_dataset(cell)
    @test dataset.uni_number == 97
    @test dataset.msg_type == 4
    @test dataset.hall_number == 90
    @test dataset.tensor_rank == 1
    @test dataset.n_operations == 32
    @test dataset.rotations == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
    ]
    @test dataset.translations ≈ [
        [0.0, 0.0, 0.0],
        [6.93889390e-18, 8.32667268e-17, 6.93889390e-17],
        [6.93889390e-18, 0.0, 6.93889390e-17],
        [0.0, 8.32667268e-17, 0.0],
        [0.75, 0.5, 0.75],
        [0.75, 0.5, 0.75],
        [0.75, 0.5, 0.75],
        [0.75, 0.5, 0.75],
        [0.5, 0.0, 0.5],
        [0.5, 8.32667268e-17, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 8.32667268e-17, 0.5],
        [0.25, 0.5, 0.25],
        [0.25, 0.5, 0.25],
        [0.25, 0.5, 0.25],
        [0.25, 0.5, 0.25],
        [0.5, 0.0, 0.0],
        [0.5, 8.32667268e-17, 6.93889390e-17],
        [0.5, 0.0, 6.93889390e-17],
        [0.5, 8.32667268e-17, 0.0],
        [0.25, 0.5, 0.75],
        [0.25, 0.5, 0.75],
        [0.25, 0.5, 0.75],
        [0.25, 0.5, 0.75],
        [0.0, 0.0, 0.5],
        [6.93889390e-18, 8.32667268e-17, 0.5],
        [6.93889390e-18, 0.0, 0.5],
        [0.0, 8.32667268e-17, 0.5],
        [0.75, 0.5, 0.25],
        [0.75, 0.5, 0.25],
        [0.75, 0.5, 0.25],
        [0.75, 0.5, 0.25],
    ]
    @test dataset.time_reversals == [
        false,
        false,
        true,
        true,
        true,
        true,
        false,
        false,
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
        false,
        false,
        true,
        true,
        true,
        true,
        false,
        false,
        false,
        false,
        true,
        true,
    ]
    @test dataset.n_atoms == 24
    @test dataset.equivalent_atoms ==
        [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 16, 16, 16, 16, 16, 16, 16
    ] .+ 1
    @test dataset.transformation_matrix == [
        -2.0 0.0 0.0
        0.0 -1.0 0.0
        1.0 0.0 1.0
    ]  # TODO: Check this
    @test dataset.origin_shift ≈ [-5.55111512e-17, -0.5, 0.0]
    @test dataset.n_std_atoms == 12
    @test dataset.std_lattice ==
        Lattice([[-3.75455, 0.0, 11.2715], [0.0, -5.7319, 0.0], [0.0, 0.0, 22.543]])
    @test dataset.std_types == [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1] .+ 1
    @test dataset.std_positions ≈ [
        [-4.44089210e-16, 0.0, 0.35125],
        [0.5, 0.5, 0.35125],
        [-2.77555756e-16, 0.0, 0.85125],
        [0.5, 0.5, 0.85125],
        [0.5, 0.5, 0.14875],
        [-2.22044605e-16, 0.0, 0.14875],
        [0.5, 0.5, 0.64875],
        [-2.22044605e-16, 0.0, 0.64875],
        [-1.11455983e-16, 0.5, -2.81892565e-18],
        [0.5, 0.0, 8.12569238e-17],
        [-2.22044605e-16, 0.5, 0.5],
        [0.5, 0.0, 0.5],
    ]
    @test dataset.std_tensors == [
        [0.12, 0.0, 9.06],
        [0.12, 0.0, 9.06],
        [-0.12, 0.0, -9.06],
        [-0.12, 0.0, -9.06],
        [-0.12, 0.0, -9.06],
        [-0.12, 0.0, -9.06],
        [0.12, 0.0, 9.06],
        [0.12, 0.0, 9.06],
        [0.58, 0.0, -1.26],
        [0.58, 0.0, -1.26],
        [-0.58, 0.0, 1.26],
        [-0.58, 0.0, 1.26],
    ]
    @test dataset.std_rotation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.primitive_lattice ≈ Lattice([
        [8.88178420e-16, 5.7319, -3.55271368e-15],
        [1.877275, -2.86595, -5.63575],
        [-7.5091, 0, 0],
    ])
end

# From https://github.com/spglib/spglib/blob/v2.6.0/test/functional/python/test_magnetic_dataset.py#L686-L718
@testset "Test primitive lattice" begin
    lattice = [[10.6949, 0, 0], [0, 6.2875, 0], [0, 0, 5.056]]
    positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.5], [0.5, 0.5, 0.5], [0.0, 0.5, 0.0]]
    atoms = [0, 0, 0, 0]
    magmoms = [[3.0, 0.4, 0.0], [-3.0, -0.4, 0.0], [-3.0, 0.4, 0.0], [3.0, -0.4, 0.0]]
    cell = SpglibCell(lattice, positions, atoms, magmoms)
    @test get_magnetic_dataset(cell).primitive_lattice ≈
        Lattice([[0, 0, -5.056], [0, -6.2875, 0], [-10.6949, 0, 0]])
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
    @testset "Test `get_symmetry_with_site_tensors`" begin
        rotations, translations, spin_flips = get_symmetry_with_site_tensors(cell, 1e-5)
        @test rotations == dataset.rotations
        @test translations == dataset.translations
        @test dataset.time_reversals == (1 .- spin_flips) ./ 2  # See https://github.com/spglib/spglib/blob/v2.1.0/test/functional/c/test_magnetic_symmetry.cpp#L156-L160
        @test sum(dataset.time_reversals) == 8  # See https://github.com/spglib/spglib/blob/f6abb97/test/functional/fortran/test_fortran_spg_get_symmetry_with_site_tensors.F90#L96
    end
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
    @test dataset.n_operations == 32  # FIXME: 96?
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
    @testset "Test `get_symmetry_with_site_tensors`" begin
        rotations, translations, spin_flips = get_symmetry_with_site_tensors(cell, 1e-5)
        @test rotations == dataset.rotations
        @test translations == dataset.translations
        @test dataset.time_reversals == (1 .- spin_flips) ./ 2  # See https://github.com/spglib/spglib/blob/v2.1.0/test/functional/c/test_magnetic_symmetry.cpp#L156-L160
        @test sum(dataset.time_reversals) == 16  # See https://github.com/spglib/spglib/blob/f6abb97/test/functional/fortran/test_fortran_spg_get_symmetry_with_site_tensors.F90#L145
    end
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
    @testset "Test `get_symmetry_with_site_tensors`" begin
        rotations, translations, spin_flips = get_symmetry_with_site_tensors(cell, 1e-5)
        @test rotations == dataset.rotations
        @test translations == dataset.translations
        @test dataset.time_reversals == (1 .- spin_flips) ./ 2  # See https://github.com/spglib/spglib/blob/v2.1.0/test/functional/c/test_magnetic_symmetry.cpp#L156-L160
    end
end

@testset "Test `get_magnetic_symmetry_from_database`" begin
    # See https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_spglib.py#L297-L300
    @testset "Test hexagonal axes: hall_number 433" begin
        rotations, translations, time_reversals = get_magnetic_symmetry_from_database(1242)
        @test rotations == [
            [
                1 0 0
                0 1 0
                0 0 1
            ],
            [
                1 0 0
                0 1 0
                0 0 1
            ],
            [
                -1 1 0
                -1 0 0
                0 0 1
            ],
            [
                0 -1 0
                1 -1 0
                0 0 1
            ],
            [
                -1 1 0
                -1 0 0
                0 0 1
            ],
            [
                -1 1 0
                -1 0 0
                0 0 1
            ],
            [
                1 0 0
                0 1 0
                0 0 1
            ],
            [
                1 0 0
                0 1 0
                0 0 1
            ],
            [
                0 -1 0
                1 -1 0
                0 0 1
            ],
            [
                -1 1 0
                -1 0 0
                0 0 1
            ],
            [
                -1 1 0
                -1 0 0
                0 0 1
            ],
            [
                0 -1 0
                1 -1 0
                0 0 1
            ],
            [
                0 -1 0
                1 -1 0
                0 0 1
            ],
            [
                -1 1 0
                -1 0 0
                0 0 1
            ],
            [
                1 0 0
                0 1 0
                0 0 1
            ],
            [
                1 0 0
                0 1 0
                0 0 1
            ],
            [
                0 -1 0
                1 -1 0
                0 0 1
            ],
            [
                0 -1 0
                1 -1 0
                0 0 1
            ],
        ]
        @test translations ≈ [
            [0.0, 0.0, 0.0],
            [0.66666667, 0.33333333, 0.33333333],
            [0.33333333, 0.66666667, 0.66666667],
            [0.0, 0.0, 0.5],
            [0.33333333, 0.66666667, 0.16666667],
            [0.66666667, 0.33333333, 0.83333333],
            [0.33333333, 0.66666667, 0.66666667],
            [0.66666667, 0.33333333, 0.83333333],
            [0.33333333, 0.66666667, 0.16666667],
            [0.66666667, 0.33333333, 0.33333333],
            [0.0, 0.0, 0.0],
            [0.66666667, 0.33333333, 0.33333333],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, 0.5],
            [0.33333333, 0.66666667, 0.16666667],
            [0.66666667, 0.33333333, 0.83333333],
            [0.33333333, 0.66666667, 0.66666667],
        ]
        @test time_reversals == [
            false,
            false,
            false,
            true,
            true,
            true,
            false,
            true,
            true,
            false,
            false,
            false,
            false,
            true,
            true,
            true,
            true,
            false,
        ]
    end
    # See https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_spglib.py#L302-L363
    @testset "Test rhombohedral axes: hall_number 434" begin
        rotations, translations, time_reversals = get_magnetic_symmetry_from_database(
            1242, 434
        )
        @test rotations == [
            [
                1 0 0
                0 1 0
                0 0 1
            ],  # x,y,z
            [
                0 0 1
                1 0 0
                0 1 0
            ],  # y,z,x
            [
                0 0 1
                1 0 0
                0 1 0
            ],  # y+1/2,z+1/2,x+1/2
            [
                0 1 0
                0 0 1
                1 0 0
            ],  # z,x,y
            [
                1 0 0
                0 1 0
                0 0 1
            ],  # x+1/2,y+1/2,z+1/2
            [
                0 1 0
                0 0 1
                1 0 0
            ],  # z+1/2,x+1/2,y+1/2
        ]
        @test translations == [
            [0, 0, 0],
            [0, 0, 0],
            [0.5, 0.5, 0.5],
            [0, 0, 0],
            [0.5, 0.5, 0.5],
            [0.5, 0.5, 0.5],
        ]
        @test time_reversals == [0, 0, 1, 0, 1, 1]
    end
end

# From https://github.com/spglib/spglib/blob/v2.1.0/test/functional/python/test_spglib.py#L257-L292
@testset "Test `get_magnetic_spacegroup_type`" begin
    @test get_magnetic_spacegroup_type(1279) ==
        MagneticSpacegroupType(1279, 1279, "156.49", "156.1.1279", 156, 1)
    @test get_magnetic_spacegroup_type(452) ==
        MagneticSpacegroupType(452, 442, "55.354", "55.2.442", 55, 2)
    @test get_magnetic_spacegroup_type(1262) ==
        MagneticSpacegroupType(1262, 1270, "151.32", "153.4.1270", 151, 4)
end
