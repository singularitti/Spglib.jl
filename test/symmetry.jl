@testset "Test `get_spacegroup_type`" begin
    # Adapted from https://github.com/unkcpz/LibSymspg.jl/blob/53d2f6d/test/test_api.jl#L7-L12
    @test get_spacegroup_type(101) == SpacegroupType(
        15,
        "C2/c",
        "I 1 1 2/a",
        "C 2/c = I 1 1 2/a",
        "C2h^6",
        101,
        "-I 2a",
        "-c3",
        "2/m",
        "C2h",
        8,
        "2/mC",
    )
    @test get_spacegroup_type(419) == SpacegroupType(
        136,
        "P4_2/mnm",
        "P 4_2/m 2_1/n 2/m",
        "P 4_2/m n m",
        "D4h^14",
        419,
        "-P 4n 2n",
        "",
        "4/mmm",
        "D4h",
        36,
        "4/mmmP",
    )
    @test get_spacegroup_type(1) ==
        SpacegroupType(1, "P1", "P 1", "P 1", "C1^1", 1, "P 1", "", "1", "C1", 1, "1P")
    @test get_spacegroup_type(525) == SpacegroupType(
        227,
        "Fd-3m",
        "F 4_1/d -3 2/m",
        "F d -3 m",
        "Oh^7",
        525,
        "F 4d 2 3 -1d",
        "1",
        "m-3m",
        "Oh",
        72,
        "m-3mF",
    )
    @test get_spacegroup_type(485) == SpacegroupType(
        191,
        "P6/mmm",
        "P 6/m 2/m 2/m",
        "P 6/m m m",
        "D6h^1",
        485,
        "-P 6 2",
        "",
        "6/mmm",
        "D6h",
        58,
        "6/mmm",
    )
end

# From https://spglib.github.io/spglib/definition.html#computing-rigid-rotation-introduced-by-idealization
@testset "Computing rigid rotation introduced by idealization" begin
    lattice = [
        [5.0759761474456697, 5.0759761474456697, 0],
        [-2.8280307701821314, 2.8280307701821314, 0],
        [0, 0, 8.57154746],
    ]
    positions = [
        [0.0, 0.84688439, 0.1203133],
        [0.0, 0.65311561, 0.6203133],
        [0.0, 0.34688439, 0.3796867],
        [0.0, 0.15311561, 0.8796867],
        [0.5, 0.34688439, 0.1203133],
        [0.5, 0.15311561, 0.6203133],
        [0.5, 0.84688439, 0.3796867],
        [0.5, 0.65311561, 0.8796867],
    ]
    atoms = fill(35, length(positions))
    cell = Cell(lattice, positions, atoms)
    dataset = get_dataset(cell, 1e-5)
    # Compared with documented results
    @test dataset.spacegroup_number == 64
    @test dataset.international_symbol == "Cmce"
    @test get_international(cell, 1e-5) == "Cmce"
    @test dataset.hall_number == 304  # Compared with Python results
    @test dataset.hall_symbol == "-C 2bc 2"  # Compared with Python results
    @test dataset.transformation_matrix ≈ [  # Compared with documented results
        1 0 0
        0 1 0
        0 0 1
    ]
    std_lattice_before_idealization =
        convert(Matrix{Float64}, Lattice(cell)) * inv(dataset.transformation_matrix)
    @testset "Test the transformation between an arbitrary system and a standardized system" begin
        @test std_lattice_before_idealization ≈ Lattice(
            [
                5.07597615 -2.82803077 0.0
                5.07597615 2.82803077 0.0
                0.0 0.0 8.57154746
            ]
        )  # Compared with Python results, the Python version is a transposed version of this
        @test std_lattice_before_idealization * dataset.transformation_matrix ≈
            Lattice(cell)
    end
    @test isempty(dataset.choice)
    @test isapprox(dataset.origin_shift, [5.55111512e-17, 0, 0]; atol=1e-16)
    @test dataset.pointgroup_symbol == "mmm"
    @test dataset.std_lattice ≈
        Lattice([[7.17851431 0.0 0.0], [0.0 3.99943947 0.0], [0.0 0.0 8.57154746]])  # Compared with Python results, the Python version is a transposed version of this
    @test dataset.std_positions ≈ [
        [0.0, 0.84688439, 0.1203133],
        [0.0, 0.65311561, 0.6203133],
        [0.0, 0.34688439, 0.3796867],
        [0.0, 0.15311561, 0.8796867],
        [0.5, 0.34688439, 0.1203133],
        [0.5, 0.15311561, 0.6203133],
        [0.5, 0.84688439, 0.3796867],
        [0.5, 0.65311561, 0.8796867],
    ]  # Compared with Python results
    @test dataset.std_types == [35, 35, 35, 35, 35, 35, 35, 35] / 35
    @test dataset.std_rotation_matrix ≈ [
        0.70710678 0.70710678 0.0
        -0.70710678 0.70710678 0.0
        0.0 0.0 1.0
    ]  # Compared with Python results
    @testset "Test the rotation of idealization" begin
        @test dataset.std_rotation_matrix ≈
            dataset.std_lattice * inv(std_lattice_before_idealization)
        @test dataset.std_lattice ≈
            dataset.std_rotation_matrix * std_lattice_before_idealization
    end
    @test dataset.mapping_to_primitive .- 1 == [0, 1, 2, 3, 0, 1, 2, 3]
    @test dataset.std_mapping_to_primitive .- 1 == [0, 1, 2, 3, 0, 1, 2, 3]
    @test dataset.wyckoffs == ['f', 'f', 'f', 'f', 'f', 'f', 'f', 'f']
    @test dataset.site_symmetry_symbols ==
        ["m..", "m..", "m..", "m..", "m..", "m..", "m..", "m.."]
    @test dataset.equivalent_atoms .- 1 == [0, 0, 0, 0, 0, 0, 0, 0]
    @test dataset.crystallographic_orbits .- 1 == [0, 0, 0, 0, 0, 0, 0, 0]
    @test dataset.primitive_lattice ≈ Lattice(
        [2.82803077, -2.82803077, 0], [1.12397269, 3.95200346, 0], [0, 0, 8.57154746]
    )
    @test all(
        dataset.translations .≈ [
            [0.0, 0.0, 0.0],
            [-1.11022302e-16, 2.46519033e-32, 0.0],
            [-1.11022302e-16, 0.5, 0.5],
            [0.0, 0.5, 0.5],
            [0.0, 0.0, 0.0],
            [-1.11022302e-16, 2.46519033e-32, 0.0],
            [-1.11022302e-16, 0.5, 0.5],
            [0.0, 0.5, 0.5],
            [0.5, 0.5, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.5, 0.0, 0.5],
        ],
    )
    python_rotations = [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
    ]
    @test dataset.rotations == python_rotations
    @test get_symmetry(cell) == (dataset.rotations, dataset.translations)
    @test get_symmetry_from_database(dataset.hall_number)[1] == dataset.rotations
    @test get_symmetry_from_database(dataset.hall_number)[2] ≈ dataset.translations
    @test get_symmetry_from_database(dataset.hall_number)[2] == [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.0, 0.5, 0.5],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.0, 0.5, 0.5],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.5, 0.0, 0.5],
    ]  # Compared with Python results
    @test get_spacegroup_type_from_symmetry(
        dataset.rotations, dataset.translations, Lattice(cell), 1e-5
    ) == SpacegroupType(
        64,
        "Cmce",
        "C 2/m 2/c 2_1/e",
        "C m c e",
        "D2h^18",
        304,
        "-C 2bc 2",
        "",
        "mmm",
        "D2h",
        19,
        "mmmC",
    )  # Compared with Python results
    @test get_hall_number_from_symmetry(dataset.rotations, dataset.translations, 1e-5) ==
        dataset.hall_number
    @test get_multiplicity(cell, 1e-5) == length(dataset.translations)
    @test get_dataset_with_hall_number(cell, dataset.hall_number) == dataset
    @test get_schoenflies(cell, 1e-5) == "D2h^18"
end

# From https://github.com/spglib/spglib/blob/ddcc153/example/python_api/example_full.py#L85-L96
@testset "Test rutile structure" begin
    lattice = [
        4 0 0
        0 4 0
        0 0 3
    ]
    positions = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.3, 0.3, 0.0],
        [0.7, 0.7, 0.0],
        [0.2, 0.8, 0.5],
        [0.8, 0.2, 0.5],
    ]
    atoms = [14, 14, 8, 8, 8, 8]
    cell = Cell(lattice, positions, atoms)
    dataset = get_dataset(cell, 1e-5)
    @test dataset.spacegroup_number == 136
    @test dataset.hall_number == 419
    @test dataset.international_symbol == "P4_2/mnm"
    @test get_international(cell, 1e-5) == "P4_2/mnm"
    @test dataset.hall_symbol == "-P 4n 2n"
    @test isempty(dataset.choice)
    @test dataset.transformation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    std_lattice_before_idealization =
        convert(Matrix{Float64}, Lattice(cell)) * inv(dataset.transformation_matrix)
    @testset "Test the transformation between an arbitrary system and a standardized system" begin
        @test std_lattice_before_idealization ≈ Lattice([
            4 0 0
            0 4 0
            0 0 3
        ])  # Compared with Python results, the Python version is a transposed version of this
        @test std_lattice_before_idealization * dataset.transformation_matrix ≈
            Lattice(cell)
    end
    @test dataset.origin_shift == [0.0, 0.0, 0.0]
    @test dataset.translations == [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
    python_rotations = [
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
    ]
    @test dataset.rotations == python_rotations
    @test get_symmetry(cell) == (dataset.rotations, dataset.translations)
    @test get_symmetry_from_database(dataset.hall_number)[1] == dataset.rotations
    @test get_symmetry_from_database(dataset.hall_number)[2] == dataset.translations
    @test dataset.wyckoffs == ['a', 'a', 'f', 'f', 'f', 'f']
    @test dataset.site_symmetry_symbols == ["m.mm", "m.mm", "m.2m", "m.2m", "m.2m", "m.2m"]
    @test dataset.crystallographic_orbits .- 1 == [0, 0, 2, 2, 2, 2]
    @test dataset.equivalent_atoms .- 1 == [0, 0, 2, 2, 2, 2]
    @test dataset.primitive_lattice ==
        Lattice([[0.0, 0.0, 3.0], [4.0, 0.0, 0.0], [0.0, 4.0, 0.0]])  # Compared with Python results, the Python version is a transposed version of this
    @test dataset.mapping_to_primitive .- 1 == 0:5
    @test dataset.std_lattice ==
        Lattice([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 3.0]])  # Compared with Python results, the Python version is a transposed version of this
    @test dataset.std_positions == [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.3, 0.3, 0.0],
        [0.7, 0.7, 0.0],
        [0.2, 0.8, 0.5],
        [0.8, 0.2, 0.5],
    ]
    @test dataset.std_rotation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @testset "Test the rotation of idealization" begin
        @test dataset.std_rotation_matrix ≈
            dataset.std_lattice * inv(std_lattice_before_idealization)
        @test dataset.std_lattice ≈
            dataset.std_rotation_matrix * std_lattice_before_idealization
    end
    @test dataset.std_types == [1, 1, 2, 2, 2, 2]  # 14, 14, 8, 8, 8, 8
    @test dataset.std_mapping_to_primitive .- 1 == [0, 1, 2, 3, 4, 5]
    @test dataset.pointgroup_symbol == "4/mmm"
    @test get_spacegroup_type_from_symmetry(
        dataset.rotations, dataset.translations, Lattice(cell), 1e-5
    ) == SpacegroupType(
        136,
        "P4_2/mnm",
        "P 4_2/m 2_1/n 2/m",
        "P 4_2/m n m",
        "D4h^14",
        419,
        "-P 4n 2n",
        "",
        "4/mmm",
        "D4h",
        36,
        "4/mmmP",
    )  # Compared with Python results
    @test get_hall_number_from_symmetry(dataset.rotations, dataset.translations, 1e-5) ==
        dataset.hall_number
    @test get_multiplicity(cell, 1e-5) == length(dataset.translations)
    @test get_dataset_with_hall_number(cell, dataset.hall_number) == dataset
    @test get_schoenflies(cell, 1e-5) == "D4h^14"
end

# From https://github.com/spglib/spglib/blob/ddcc153/example/python_api/example_full.py#L98-L109
@testset "Test distorted rutile structure" begin
    lattice = [
        3.97 0.00 0.00
        0.00 4.03 0.00
        0.00 0.00 3.00
    ]
    positions = [
        [0.0, 0.0, 0.0],
        [0.5001, 0.5, 0.5],
        [0.3, 0.3, 0.0],
        [0.7, 0.7, 0.002],
        [0.2, 0.8, 0.5],
        [0.8, 0.2, 0.5],
    ]
    atoms = [14, 14, 8, 8, 8, 8]
    cell = Cell(lattice, positions, atoms)
    dataset = get_dataset(cell, 1e-5)
    # These results are compared with Python's spglib results.
    @test dataset.spacegroup_number == 1
    @test dataset.hall_number == 1
    @test dataset.hall_symbol == "P 1"
    @test dataset.international_symbol == "P1"
    @test get_international(cell, 1e-5) == "P1"
    @test isempty(dataset.choice)
    @test dataset.transformation_matrix == [
        0 0 1
        1 0 0
        0 1 0
    ]
    std_lattice_before_idealization =
        convert(Matrix{Float64}, Lattice(cell)) * inv(dataset.transformation_matrix)
    @testset "Test the transformation between an arbitrary system and a standardized system" begin
        @test std_lattice_before_idealization ≈
            Lattice([[0.0, 0.0, 3.0], [3.97, 0.0, 0.0], [0.0, 4.03, 0.0]])  # Compared with Python results, the Python version is a transposed version of this
        @test std_lattice_before_idealization * dataset.transformation_matrix ≈
            Lattice(cell)
    end
    @test dataset.origin_shift == [0.0, 0.0, 0.0]
    @test only(dataset.rotations) == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.translations == [[0.0, 0.0, 0.0]]
    @test dataset.wyckoffs == ['a', 'a', 'a', 'a', 'a', 'a']
    @test dataset.site_symmetry_symbols == ["1", "1", "1", "1", "1", "1"]
    @test dataset.crystallographic_orbits .- 1 == [0, 1, 2, 3, 4, 5]
    @test dataset.equivalent_atoms .- 1 == [0, 1, 2, 3, 4, 5]
    @test dataset.primitive_lattice ==
        Lattice([[0.0, 0.0, 3.0], [3.97, 0.0, 0.0], [0.0, 4.03, 0.0]])
    @test dataset.mapping_to_primitive .- 1 == [0, 1, 2, 3, 4, 5]
    @test dataset.std_lattice ≈ Lattice([
        [3, 0, 0], [2.4309239e-16, 3.97, 0], [2.4676633e-16, 2.4676633e-16, 4.03]
    ])
    @test dataset.std_positions == [
        [0.0, 0.0, 0.0],
        [0.5, 0.5001, 0.5],
        [0.0, 0.3, 0.3],
        [0.002, 0.7, 0.7],
        [0.5, 0.2, 0.8],
        [0.5, 0.8, 0.2],
    ]
    @test dataset.std_types == [1, 1, 2, 2, 2, 2]
    @test dataset.std_rotation_matrix == [
        0 0 1
        1 0 0
        0 1 0
    ]
    @testset "Test the rotation of idealization" begin
        @test dataset.std_rotation_matrix ≈
            dataset.std_lattice * inv(std_lattice_before_idealization)
        @test dataset.std_rotation_matrix ≈ [
            6.123234e-17 6.123234e-17 1
            1 6.123234e-17 0
            0 1 0
        ]  # Python result
        @test dataset.std_lattice ≈
            dataset.std_rotation_matrix * std_lattice_before_idealization
    end
    @test dataset.std_mapping_to_primitive .- 1 == [0, 1, 2, 3, 4, 5]
    @test dataset.pointgroup_symbol == "1"
    @test get_symmetry(cell) == (dataset.rotations, dataset.translations)
    @test get_symmetry_from_database(dataset.hall_number)[1] == dataset.rotations
    @test get_symmetry_from_database(dataset.hall_number)[2] == dataset.translations
    @test get_spacegroup_type_from_symmetry(
        dataset.rotations, dataset.translations, Lattice(cell), 1e-5
    ) == SpacegroupType(1, "P1", "P 1", "P 1", "C1^1", 1, "P 1", "", "1", "C1", 1, "1P")  # Compared with Python results
    @test get_hall_number_from_symmetry(dataset.rotations, dataset.translations, 1e-5) ==
        dataset.hall_number
    @test get_multiplicity(cell, 1e-5) == length(dataset.translations)
    @test get_dataset_with_hall_number(cell, dataset.hall_number) == dataset
    @test get_schoenflies(cell, 1e-5) == "C1^1"
end

# From https://github.com/unkcpz/LibSymspg.jl/blob/53d2f6d/test/test_api.jl#L14-L32
@testset "Test on P-3m1 space group" begin
    lattice = [
        4.0 2.0 0.0
        0.0 3.4641 0.0
        0.0 0.0 12.0
    ]
    positions = [[0.0, 0.0, 0.0], [1 / 3, 1 / 3, 1 / 3]]
    atoms = [1, 1]
    cell = Cell(lattice, positions, atoms)
    dataset = get_dataset(cell, 1e-3)
    # Compared results with Python spglib
    @test dataset.spacegroup_number == 164
    @test dataset.hall_number == 456
    @test dataset.hall_symbol == "-P 3 2\""
    @test dataset.international_symbol == "P-3m1"
    @test get_international(cell, 1e-3) == "P-3m1"
    @test isempty(dataset.choice)
    @test dataset.transformation_matrix ≈ [
        1 0 0
        1 1 0
        0 0 1
    ]
    std_lattice_before_idealization =
        convert(Matrix{Float64}, Lattice(cell)) * inv(dataset.transformation_matrix)
    @testset "Test the transformation between an arbitrary system and a standardized system" begin
        @test std_lattice_before_idealization ≈ Lattice([
            2.0 2.0 0.0
            -3.4641 3.4641 0.0
            0.0 0.0 12.0
        ])  # Compared with Python results, the Python version is a transposed version of this
        @test std_lattice_before_idealization * dataset.transformation_matrix ≈
            Lattice(cell)
    end
    @test dataset.origin_shift ≈ [1 / 3, 2 / 3, 1 / 3]
    @test dataset.n_operations == 12
    python_rotations = [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 -1 0; 1 0 0; 0 0 1],
        [1 1 0; -1 0 0; 0 0 -1],
        [0 1 0; -1 -1 0; 0 0 1],
        [0 -1 0; 1 1 0; 0 0 -1],
        [1 1 0; 0 -1 0; 0 0 -1],
        [-1 -1 0; 0 1 0; 0 0 1],
        [0 -1 0; -1 0 0; 0 0 -1],
        [0 1 0; 1 0 0; 0 0 1],
        [-1 0 0; 1 1 0; 0 0 -1],
        [1 0 0; -1 -1 0; 0 0 1],
    ]
    @test dataset.rotations == python_rotations
    @test all(
        dataset.translations .≈ [
            [0.0, 0.0, 0.0],
            [1 / 3, 1 / 3, 1 / 3],
            [0.0, 0.0, 0.0],
            [1 / 3, 1 / 3, 1 / 3],
            [0.0, 0.0, 0.0],
            [1 / 3, 1 / 3, 1 / 3],
            [1 / 3, 1 / 3, 1 / 3],
            [0.0, 0.0, 0.0],
            [1 / 3, 1 / 3, 1 / 3],
            [0.0, 0.0, 0.0],
            [1 / 3, 1 / 3, 1 / 3],
            [0.0, 0.0, 0.0],
        ],
    )
    @test size(dataset.rotations) == size(dataset.translations) == (12,)
    @test get_symmetry(cell) == (dataset.rotations, dataset.translations)
    @test get_symmetry_from_database(dataset.hall_number)[1] == [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [0 -1 0; 1 -1 0; 0 0 1],
        [0 1 0; -1 1 0; 0 0 -1],
        [-1 1 0; -1 0 0; 0 0 1],
        [1 -1 0; 1 0 0; 0 0 -1],
        [0 1 0; 1 0 0; 0 0 -1],
        [0 -1 0; -1 0 0; 0 0 1],
        [1 -1 0; 0 -1 0; 0 0 -1],
        [-1 1 0; 0 1 0; 0 0 1],
        [-1 0 0; -1 1 0; 0 0 -1],
        [1 0 0; 1 -1 0; 0 0 1],
    ]  # Compared with Python results
    @test get_symmetry_from_database(dataset.hall_number)[2] == [
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
    @test dataset.wyckoffs == ['d', 'd']
    @test dataset.site_symmetry_symbols == ["3m.", "3m."]
    @test dataset.equivalent_atoms .- 1 == [0, 0]
    @test dataset.crystallographic_orbits .- 1 == [0, 0]
    @test dataset.primitive_lattice ==
        Lattice([[2.0, -3.4641, 0], [-2.0, -3.4641, 0], [0.0, 0.0, -12.0]])
    @test dataset.mapping_to_primitive .- 1 == [0, 1]
    @test dataset.std_lattice ≈
        Lattice([[3.9999986, 0, 0], [-1.9999993, 3.4641004, 0], [0, 0, 12.0]])
    @test dataset.std_types == [1, 1]
    @test dataset.std_positions ≈ [[1 / 3, 2 / 3, 1 / 3], [2 / 3, 1 / 3, 2 / 3]]
    @test dataset.std_rotation_matrix ≈ [
        0.50000017 -0.8660253 0
        0.8660253 0.50000017 0
        0 0 1
    ]
    @testset "Test the rotation of idealization" begin
        @test isapprox(
            dataset.std_rotation_matrix,
            dataset.std_lattice * inv(std_lattice_before_idealization);
            atol=1e-6,
        )
        @test isapprox(
            dataset.std_lattice,
            dataset.std_rotation_matrix * std_lattice_before_idealization;
            atol=5e-6,
        )
    end
    @test dataset.std_mapping_to_primitive .- 1 == [0, 1]
    @test dataset.pointgroup_symbol == "-3m"
    @test get_spacegroup_type_from_symmetry(
        dataset.rotations, dataset.translations, Lattice(cell), 1e-5
    ) == SpacegroupType(
        164,
        "P-3m1",
        "P -3 2/m 1",
        "P -3 m 1",
        "D3d^3",
        456,
        "-P 3 2\"",
        "",
        "-3m",
        "D3d",
        49,
        "-3m1P",
    )  # Compared with Python results
    @test get_hall_number_from_symmetry(dataset.rotations, dataset.translations, 1e-5) ==
        dataset.hall_number
    @test get_multiplicity(cell, 1e-5) == length(dataset.translations)
    @test get_dataset_with_hall_number(cell, dataset.hall_number) == dataset
    @test get_schoenflies(cell, 1e-5) == "D3d^3"
end

# From https://github.com/spglib/spglib/blob/ddcc153/example/python_api/example_full.py#L43-L59
@testset "Test silicon structure" begin
    lattice = Lattice([[4, 0, 0], [0, 4, 0], [0, 0, 4]])
    positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ]
    atoms = [14, 14, 14, 14, 14, 14, 14, 14]
    cell = Cell(lattice, positions, atoms)
    dataset = get_dataset(cell, 1e-5)
    # Compared with Python results
    @test dataset.spacegroup_number == 227
    @test dataset.hall_number == 525
    @test dataset.international_symbol == "Fd-3m"
    @test dataset.hall_symbol == "F 4d 2 3 -1d"
    @test get_international(cell, 1e-5) == "Fd-3m"
    @test dataset.choice == "1"
    @test dataset.transformation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    std_lattice_before_idealization =
        convert(Matrix{Float64}, Lattice(cell)) * inv(dataset.transformation_matrix)
    @testset "Test the transformation between an arbitrary system and a standardized system" begin
        @test std_lattice_before_idealization ≈ Lattice([
            4 0 0
            0 4 0
            0 0 4
        ])  # Compared with Python results, the Python version is a transposed version of this
        @test std_lattice_before_idealization * dataset.transformation_matrix ≈
            Lattice(cell)
    end
    @test dataset.origin_shift == [0, 0, 0]
    @test size(dataset.rotations) == (192,)
    @test dataset.translations == [
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.75, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.25, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.75, 0.75],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.0, 0.0],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
        [0.5, 0.0, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.0, 0.0],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.0, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.75, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.0, 0.5],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.75, 0.25, 0.75],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.0, 0.0, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.25, 0.75],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.75],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.0, 0.5],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.0, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
    ]  # Compared with Python results
    @test get_symmetry(cell) == (dataset.rotations, dataset.translations)
    @test get_symmetry_from_database(dataset.hall_number)[1] == dataset.rotations
    @test get_symmetry_from_database(dataset.hall_number)[2] == [
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.75, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.25, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.75, 0.75],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.0, 0.0],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.0, 0.0, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.25, 0.75],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.75],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.0, 0.5],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.0, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.75, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.0, 0.5],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.75, 0.25, 0.75],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
        [0.5, 0.0, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
        [0.5, 0.0, 0.5],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.0, 0.0],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
        [0.0, 0.5, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.0, 0.0, 0.0],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.0, 0.0],
        [0.75, 0.75, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
    ]  # Compared with Python results
    @test dataset.wyckoffs == ['a', 'a', 'a', 'a', 'a', 'a', 'a', 'a']
    @test dataset.site_symmetry_symbols ==
        ["-43m", "-43m", "-43m", "-43m", "-43m", "-43m", "-43m", "-43m"]
    @test dataset.equivalent_atoms .- 1 == [0, 0, 0, 0, 0, 0, 0, 0]
    @test dataset.crystallographic_orbits .- 1 == [0, 0, 0, 0, 0, 0, 0, 0]
    @test dataset.primitive_lattice ==
        Lattice([[2.0, -2.0, 0.0], [-2.0, -0.0, -2.0], [2.0, 2.0, 0.0]])
    @test dataset.mapping_to_primitive .- 1 == [0, 0, 0, 0, 1, 1, 1, 1]
    @test dataset.std_lattice ==
        Lattice([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]])
    @test dataset.std_types == [14, 14, 14, 14, 14, 14, 14, 14] / 14
    @test dataset.std_positions == [
        [0.0, 0.0, 0.0],
        [0.25, 0.75, 0.75],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.5, 0.0, 0.5],
        [0.75, 0.75, 0.25],
        [0.5, 0.5, 0.0],
        [0.75, 0.25, 0.75],
    ]
    @test dataset.std_rotation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @testset "Test the rotation of idealization" begin
        @test dataset.std_rotation_matrix ≈
            dataset.std_lattice * inv(std_lattice_before_idealization)
        @test dataset.std_lattice ≈
            dataset.std_rotation_matrix * std_lattice_before_idealization
    end
    @test dataset.std_mapping_to_primitive .- 1 == [0, 1, 0, 1, 0, 1, 0, 1]
    @test dataset.pointgroup_symbol == "m-3m"
    @test get_spacegroup_type_from_symmetry(
        dataset.rotations, dataset.translations, Lattice(cell), 1e-5
    ) == SpacegroupType(
        227,
        "Fd-3m",
        "F 4_1/d -3 2/m",
        "F d -3 m",
        "Oh^7",
        525,
        "F 4d 2 3 -1d",
        "1",
        "m-3m",
        "Oh",
        72,
        "m-3mF",
    )  # Compared with Python results
    @test get_hall_number_from_symmetry(dataset.rotations, dataset.translations, 1e-5) ==
        dataset.hall_number
    @test get_multiplicity(cell, 1e-5) == length(dataset.translations)
    @test get_dataset_with_hall_number(cell, dataset.hall_number) == dataset
    @test get_schoenflies(cell, 1e-5) == "Oh^7"
end

# From https://github.com/spglib/spglib/blob/ddcc153/example/python_api/example_full.py#L79-L83
@testset "Test primitive silicon structure" begin
    lattice = [[0, 2, 2], [2, 0, 2], [2, 2, 0]]
    positions = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]
    atoms = [14, 14]
    cell = Cell(lattice, positions, atoms)
    dataset = get_dataset(cell, 1e-5)
    # Compared with Python results
    @test dataset.spacegroup_number == 227
    @test dataset.hall_number == 525
    @test dataset.international_symbol == "Fd-3m"
    @test dataset.hall_symbol == "F 4d 2 3 -1d"
    @test get_international(cell, 1e-5) == "Fd-3m"
    @test dataset.choice == "1"
    @test dataset.transformation_matrix == [
        0.0 0.5 0.5
        0.5 0.0 0.5
        0.5 0.5 0.0
    ]
    std_lattice_before_idealization =
        convert(Matrix{Float64}, Lattice(cell)) * inv(dataset.transformation_matrix)
    @testset "Test the transformation between an arbitrary system and a standardized system" begin
        @test std_lattice_before_idealization ≈
            Lattice([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]])  # Compared with Python results, the Python version is a transposed version of this
        @test std_lattice_before_idealization * dataset.transformation_matrix ≈
            Lattice(cell)
    end
    @test dataset.origin_shift == [0, 0, 1 / 2]
    python_rotations = [
        [1 0 0; 0 1 0; 0 0 1],
        [1 1 1; 0 0 -1; -1 0 0],
        [0 1 0; 1 0 0; -1 -1 -1],
        [0 0 -1; 1 1 1; 0 -1 0],
        [-1 -1 -1; 0 0 1; 0 1 0],
        [0 -1 0; -1 0 0; 0 0 -1],
        [0 0 1; -1 -1 -1; 1 0 0],
        [-1 0 0; 0 -1 0; 1 1 1],
        [0 0 1; 1 0 0; 0 1 0],
        [-1 0 0; 1 1 1; 0 0 -1],
        [-1 -1 -1; 0 1 0; 1 0 0],
        [0 -1 0; 0 0 -1; 1 1 1],
        [0 1 0; -1 -1 -1; 0 0 1],
        [0 0 -1; 0 -1 0; -1 0 0],
        [1 0 0; 0 0 1; -1 -1 -1],
        [1 1 1; -1 0 0; 0 -1 0],
        [0 1 0; 0 0 1; 1 0 0],
        [0 0 -1; -1 0 0; 1 1 1],
        [1 0 0; -1 -1 -1; 0 1 0],
        [1 1 1; 0 -1 0; 0 0 -1],
        [0 0 1; 0 1 0; -1 -1 -1],
        [-1 0 0; 0 0 -1; 0 -1 0],
        [-1 -1 -1; 1 0 0; 0 0 1],
        [0 -1 0; 1 1 1; -1 0 0],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 -1 -1; 0 0 1; 1 0 0],
        [0 -1 0; -1 0 0; 1 1 1],
        [0 0 1; -1 -1 -1; 0 1 0],
        [1 1 1; 0 0 -1; 0 -1 0],
        [0 1 0; 1 0 0; 0 0 1],
        [0 0 -1; 1 1 1; -1 0 0],
        [1 0 0; 0 1 0; -1 -1 -1],
        [0 0 -1; -1 0 0; 0 -1 0],
        [1 0 0; -1 -1 -1; 0 0 1],
        [1 1 1; 0 -1 0; -1 0 0],
        [0 1 0; 0 0 1; -1 -1 -1],
        [0 -1 0; 1 1 1; 0 0 -1],
        [0 0 1; 0 1 0; 1 0 0],
        [-1 0 0; 0 0 -1; 1 1 1],
        [-1 -1 -1; 1 0 0; 0 1 0],
        [0 -1 0; 0 0 -1; -1 0 0],
        [0 0 1; 1 0 0; -1 -1 -1],
        [-1 0 0; 1 1 1; 0 -1 0],
        [-1 -1 -1; 0 1 0; 0 0 1],
        [0 0 -1; 0 -1 0; 1 1 1],
        [1 0 0; 0 0 1; 0 1 0],
        [1 1 1; -1 0 0; 0 0 -1],
        [0 1 0; -1 -1 -1; 1 0 0],
    ]  # Compared with Python results
    @test dataset.rotations == python_rotations
    @test dataset.translations == [
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.0, 0.0],
    ]
    @test get_symmetry(cell) == (dataset.rotations, dataset.translations)
    @test size(get_symmetry_from_database(dataset.hall_number)[1]) ==
        size(get_symmetry_from_database(dataset.hall_number)[2]) ==
        (192,)  # Symmetry breaking
    @test dataset.wyckoffs == ['b', 'b']
    @test dataset.site_symmetry_symbols == ["-43m", "-43m"]
    @test dataset.crystallographic_orbits .- 1 == [0, 0]
    @test dataset.equivalent_atoms .- 1 == [0, 0]
    @test dataset.primitive_lattice ==
        Lattice([[2.0, -2.0, 0.0], [-2.0, -0.0, -2.0], [2.0, 2.0, 0.0]])
    @test dataset.mapping_to_primitive .- 1 == [0, 1]
    @test dataset.std_lattice ==
        Lattice([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]])
    @test dataset.std_types == [14, 14, 14, 14, 14, 14, 14, 14] / 14
    @test dataset.std_positions == [
        [0.0, 0.0, 0.5],
        [0.25, 0.75, 0.25],
        [0.0, 0.5, 0.0],
        [0.25, 0.25, 0.75],
        [0.5, 0.0, 0.0],
        [0.75, 0.75, 0.75],
        [0.5, 0.5, 0.5],
        [0.75, 0.25, 0.25],
    ]
    @test dataset.std_rotation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @testset "Test the rotation of idealization" begin
        @test dataset.std_rotation_matrix ≈
            dataset.std_lattice * inv(std_lattice_before_idealization)
        @test dataset.std_lattice ≈
            dataset.std_rotation_matrix * std_lattice_before_idealization
    end
    @test dataset.std_mapping_to_primitive .- 1 == [0, 1, 0, 1, 0, 1, 0, 1]
    @test dataset.pointgroup_symbol == "m-3m"
    @test get_spacegroup_type_from_symmetry(
        dataset.rotations, dataset.translations, Lattice(cell), 1e-5
    ) == SpacegroupType(
        227,
        "Fd-3m",
        "F 4_1/d -3 2/m",
        "F d -3 m",
        "Oh^7",
        525,
        "F 4d 2 3 -1d",
        "1",
        "m-3m",
        "Oh",
        72,
        "m-3mF",
    )  # Compared with Python results
    @test get_hall_number_from_symmetry(dataset.rotations, dataset.translations, 1e-5) ==
        dataset.hall_number
    @test get_multiplicity(cell, 1e-5) == length(dataset.translations)
    @test get_dataset_with_hall_number(cell, dataset.hall_number) == dataset
    @test get_schoenflies(cell, 1e-5) == "Oh^7"
end

# From https://github.com/spglib/spglib/blob/ddcc153/example/python_api/example_full.py#L111-L117
@testset "Test MgB₂ structure" begin
    a = 3.07
    c = 3.52
    lattice = [[a, 0, 0], [-a / 2, a / 2 * sqrt(3), 0], [0, 0, c]]
    positions = [[0, 0, 0], [1 / 3, 2 / 3, 1 / 2], [2 / 3, 1 / 3, 1 / 2]]
    atoms = [12, 5, 5]
    cell = Cell(lattice, positions, atoms)
    dataset = get_dataset(cell, 1e-5)
    @test dataset.spacegroup_number == 191  # Compared with C results
    @test dataset.hall_number == 485
    @test dataset.hall_symbol == "-P 6 2"
    @test dataset.international_symbol == "P6/mmm"
    @test get_international(cell) == "P6/mmm"
    @test dataset.pointgroup_symbol == "6/mmm"
    @test isempty(dataset.choice)
    @test dataset.transformation_matrix ≈ [
        1 5.55111512e-17 0
        0 1 0
        0 0 1
    ]
    std_lattice_before_idealization =
        convert(Matrix{Float64}, Lattice(cell)) * inv(dataset.transformation_matrix)
    @testset "Test the transformation between an arbitrary system and a standardized system" begin
        @test std_lattice_before_idealization ≈
            Lattice([[3.07, 0, 0], [-1.535, 2.65869799, 0], [0, 0, 3.52]])  # Compared with Python results, the Python version is a transposed version of this
        @test std_lattice_before_idealization * dataset.transformation_matrix ≈
            Lattice(cell)
    end
    @test dataset.origin_shift == [0, 0, 0]
    python_rotations = [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [1 -1 0; 1 0 0; 0 0 1],
        [-1 1 0; -1 0 0; 0 0 -1],
        [0 -1 0; 1 -1 0; 0 0 1],
        [0 1 0; -1 1 0; 0 0 -1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
        [-1 1 0; -1 0 0; 0 0 1],
        [1 -1 0; 1 0 0; 0 0 -1],
        [0 1 0; -1 1 0; 0 0 1],
        [0 -1 0; 1 -1 0; 0 0 -1],
        [0 -1 0; -1 0 0; 0 0 -1],
        [0 1 0; 1 0 0; 0 0 1],
        [-1 0 0; -1 1 0; 0 0 -1],
        [1 0 0; 1 -1 0; 0 0 1],
        [-1 1 0; 0 1 0; 0 0 -1],
        [1 -1 0; 0 -1 0; 0 0 1],
        [0 1 0; 1 0 0; 0 0 -1],
        [0 -1 0; -1 0 0; 0 0 1],
        [1 0 0; 1 -1 0; 0 0 -1],
        [-1 0 0; -1 1 0; 0 0 1],
        [1 -1 0; 0 -1 0; 0 0 -1],
        [-1 1 0; 0 1 0; 0 0 1],
    ]
    @test dataset.translations == [
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
    ]  # Compared with Python results
    @test get_symmetry(cell, 1e-5) == (dataset.rotations, dataset.translations)
    @test get_symmetry_from_database(dataset.hall_number)[1] == dataset.rotations
    @test get_symmetry_from_database(dataset.hall_number)[2] == dataset.translations
    @test dataset.wyckoffs == ['a', 'd', 'd']
    @test dataset.site_symmetry_symbols == ["6/mmm", "-6m2", "-6m2"]
    @test dataset.crystallographic_orbits .- 1 == [0, 1, 1]
    @test dataset.equivalent_atoms .- 1 == [0, 1, 1]
    @test dataset.primitive_lattice ≈
        Lattice([[3.07, 0.0, 0.0], [-1.535, 2.65869799, 0.0], [0.0, 0.0, 3.52]])
    @test dataset.mapping_to_primitive .- 1 == 0:2
    @test dataset.std_lattice ≈
        Lattice([[3.07, 0.0, 0.0], [-1.535, 2.65869799, 0.0], [0.0, 0.0, 3.52]])
    @test dataset.std_types == [1, 2, 2]
    @test dataset.std_positions ≈
        [[0.0, 0.0, 0.0], [1 / 3, 2 / 3, 0.5], [2 / 3, 1 / 3, 0.5]]
    @test dataset.std_rotation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @testset "Test the rotation of idealization" begin
        @test dataset.std_rotation_matrix ≈
            dataset.std_lattice * inv(std_lattice_before_idealization)
        @test dataset.std_lattice ≈
            dataset.std_rotation_matrix * std_lattice_before_idealization
    end
    @test dataset.std_mapping_to_primitive .- 1 == 0:2
    @test get_spacegroup_type_from_symmetry(
        dataset.rotations, dataset.translations, Lattice(cell), 1e-5
    ) == SpacegroupType(
        191,
        "P6/mmm",
        "P 6/m 2/m 2/m",
        "P 6/m m m",
        "D6h^1",
        485,
        "-P 6 2",
        "",
        "6/mmm",
        "D6h",
        58,
        "6/mmm",
    )  # Compared with Python results
    @test get_hall_number_from_symmetry(dataset.rotations, dataset.translations, 1e-5) ==
        dataset.hall_number
    @test get_multiplicity(cell, 1e-5) == length(dataset.translations)
    @test get_dataset_with_hall_number(cell, dataset.hall_number) == dataset
    @test get_schoenflies(cell, 1e-5) == "D6h^1"
end

# Example is from here: https://github.com/spglib/spglib/blob/ddcc153/example/python_api/example.py
@testset "Test wurtzite structure (P6_3mc)" begin
    lattice = [[3.111, 0, 0], [-1.5555, 2.6942050311733885, 0], [0, 0, 4.988]]  # Note this is different from C
    positions = [
        [1.0 / 3, 2.0 / 3, 0.0],
        [2.0 / 3, 1.0 / 3, 0.5],
        [1.0 / 3, 2.0 / 3, 0.6181],
        [2.0 / 3, 1.0 / 3, 0.1181],
    ]
    atoms = [1, 1, 2, 2]
    cell = Cell(lattice, positions, atoms)
    dataset = get_dataset(cell, 1e-5)
    @test get_international(cell) == "P6_3mc"
    @test dataset.hall_number == 480
    @test dataset.pointgroup_symbol == "6mm"
    @test dataset.spacegroup_number == 186  # Compared with C results
    @test dataset.hall_symbol == "P 6c -2c"
    @test isempty(dataset.choice)
    @test dataset.transformation_matrix ≈ [
        1 -5.55111512e-17 0
        0 1 0
        0 0 1
    ]
    std_lattice_before_idealization =
        convert(Matrix{Float64}, Lattice(cell)) * inv(dataset.transformation_matrix)
    @testset "Test the transformation between an arbitrary system and a standardized system" begin
        @test std_lattice_before_idealization ≈ Lattice([
            3.111 -1.5555 0
            0 2.69420503 0
            0 0 4.988
        ])  # Compared with Python results, the Python version is a transposed version of this
        @test std_lattice_before_idealization * dataset.transformation_matrix ≈
            Lattice(cell)
    end
    @test dataset.origin_shift ≈ [-5.55111512e-17, 0, 0]  # Compared with Python results
    python_rotations = [
        [1 0 0; 0 1 0; 0 0 1],
        [1 -1 0; 1 0 0; 0 0 1],
        [0 -1 0; 1 -1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [-1 1 0; -1 0 0; 0 0 1],
        [0 1 0; -1 1 0; 0 0 1],
        [0 1 0; 1 0 0; 0 0 1],
        [1 0 0; 1 -1 0; 0 0 1],
        [1 -1 0; 0 -1 0; 0 0 1],
        [0 -1 0; -1 0 0; 0 0 1],
        [-1 0 0; -1 1 0; 0 0 1],
        [-1 1 0; 0 1 0; 0 0 1],
    ]
    @test dataset.rotations == python_rotations
    @test dataset.translations ≈ [
        [0, 0, 0],
        [3.08148791e-33, -5.55111512e-17, 0.5],
        [5.55111512e-17, -5.55111512e-17, 0],
        [1.11022302e-16, 0, 0.5],
        [1.11022302e-16, 5.55111512e-17, 0],
        [5.55111512e-17, 5.55111512e-17, 0.5],
        [5.55111512e-17, -5.55111512e-17, 0.5],
        [3.08148791e-33, -5.55111512e-17, 0],
        [0, 0, 0.5],
        [5.55111512e-17, 5.55111512e-17, 0],
        [1.11022302e-16, 5.55111512e-17, 0.5],
        [1.11022302e-16, 0, 0],
    ]  # Compared with Python results
    @test get_symmetry(cell) == (dataset.rotations, dataset.translations)
    @test get_symmetry_from_database(dataset.hall_number)[1] == dataset.rotations
    @test get_symmetry_from_database(dataset.hall_number)[2] ≈ dataset.translations
    @test get_symmetry_from_database(dataset.hall_number)[2] == [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5],
        [0.0, 0.0, 0.5],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5],
        [0.0, 0.0, 0.0],
    ]  # Compared with Python results
    @test dataset.wyckoffs == ['b', 'b', 'b', 'b']
    @test dataset.site_symmetry_symbols == ["3m.", "3m.", "3m.", "3m."]
    @test dataset.crystallographic_orbits .- 1 == [0, 0, 2, 2]
    @test dataset.equivalent_atoms .- 1 == [0, 0, 2, 2]
    @test dataset.primitive_lattice ≈
        Lattice([[3.111, 0, 0], [-1.5555, 2.69420503, 0], [0, 0, 4.988]])  # Compared with Python results, the Python version is a transposed version of this
    @test dataset.mapping_to_primitive .- 1 == 0:3  # FIXME: should I +1?
    @test dataset.std_lattice ≈
        Lattice([[3.111, 0, 0], [-1.5555, 2.69420503, 0], [0, 0, 4.988]])  # Compared with Python results, the Python version is a transposed version of this
    @test dataset.std_positions ≈ [
        [0.33333333, 0.66666667, 0],
        [0.66666667, 0.33333333, 0.5],
        [0.33333333, 0.66666667, 0.6181],
        [0.66666667, 0.33333333, 0.1181],
    ]  # Compared with Python results
    @test dataset.std_rotation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @testset "Test the rotation of idealization" begin
        @test dataset.std_rotation_matrix ≈
            dataset.std_lattice * inv(std_lattice_before_idealization)
        @test dataset.std_lattice ≈
            dataset.std_rotation_matrix * std_lattice_before_idealization
    end
    @test dataset.std_mapping_to_primitive .- 1 == 0:3  # FIXME: should I +1?
    @test get_spacegroup_type_from_symmetry(
        dataset.rotations, dataset.translations, Lattice(cell), 1e-5
    ) == SpacegroupType(
        186,
        "P6_3mc",
        "P 6_3 m c",
        "P 6_3 m c",
        "C6v^4",
        480,
        "P 6c -2c",
        "",
        "6mm",
        "C6v",
        55,
        "6mmP",
    )  # Compared with Python results
    @test get_hall_number_from_symmetry(dataset.rotations, dataset.translations, 1e-5) ==
        dataset.hall_number
    @test get_multiplicity(cell, 1e-5) == length(dataset.translations)
    @test get_dataset_with_hall_number(cell, dataset.hall_number) == dataset
    @test get_schoenflies(cell, 1e-5) == "C6v^4"
end
