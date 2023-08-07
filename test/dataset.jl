# From https://spglib.github.io/spglib/definition.html#computing-rigid-rotation-introduced-by-idealization
@testset "Computing rigid rotation introduced by idealization" begin
    lattice = [
        [5.0759761474456697, 5.0759761474456697, 0],  # a
        [-2.8280307701821314, 2.8280307701821314, 0],  # b
        [0, 0, 8.57154746],  # c
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
    numbers = fill(35, length(positions))
    cell = Cell(lattice, positions, numbers)
    dataset = get_dataset(cell, 1e-5)
    # Compared with documented results
    @test get_spacegroup_number(cell) == 64
    @test dataset.spacegroup_number == 64
    @test dataset.international_symbol == "Cmce"
    @test get_international(cell, 1e-5) == "Cmce"
    @test dataset.hall_number == 304  # Compared with Python results
    @test dataset.hall_symbol == "-C 2bc 2"  # Compared with Python results
    @test dataset.transformation_matrix == [  # Compared with documented results
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.std_lattice ≈ Lattice([  # Compared with documented results
        7.17851431 0 0
        0 3.99943947 0
        0 0 8.57154746
    ])
    # Compared with Python results
    @test isempty(dataset.choice)
    @test isapprox(dataset.origin_shift, [5.55111512e-17, 0, 0]; atol=1e-16)
    @test dataset.pointgroup_symbol == "mmm"
    @test dataset.std_types == [35, 35, 35, 35, 35, 35, 35, 35] / 35
    @test dataset.std_positions ≈ [
        [0.0, 0.84688439, 0.1203133],
        [0.0, 0.65311561, 0.6203133],
        [0.0, 0.34688439, 0.3796867],
        [0.0, 0.15311561, 0.8796867],
        [0.5, 0.34688439, 0.1203133],
        [0.5, 0.15311561, 0.6203133],
        [0.5, 0.84688439, 0.3796867],
        [0.5, 0.65311561, 0.8796867],
    ]
    @test dataset.std_rotation_matrix ≈ [
        0.70710678 0.70710678 0.0
        -0.70710678 0.70710678 0.0
        0.0 0.0 1.0
    ]
    @test dataset.mapping_to_primitive == [0, 1, 2, 3, 0, 1, 2, 3]
    @test dataset.std_mapping_to_primitive == [0, 1, 2, 3, 0, 1, 2, 3]
    @test dataset.wyckoffs == ['f', 'f', 'f', 'f', 'f', 'f', 'f', 'f']
    @test dataset.site_symmetry_symbols ==
        ["m..", "m..", "m..", "m..", "m..", "m..", "m..", "m.."]
    @test dataset.equivalent_atoms == [0, 0, 0, 0, 0, 0, 0, 0]
    @test dataset.crystallographic_orbits == [0, 0, 0, 0, 0, 0, 0, 0]
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
    @test all(
        map(zip(dataset.rotations, python_rotations)) do (rotation, python_rotation)
            rotation == python_rotation
        end,
    )
    @test get_symmetry(cell) == (dataset.rotations, dataset.translations)
end

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
    types = [14, 14, 8, 8, 8, 8]
    rutile = Cell(lattice, positions, types)
    dataset = get_dataset(rutile, 1e-5)
    @test dataset.spacegroup_number == 136
    @test dataset.site_symmetry_symbols == ["m.mm", "m.mm", "m.2m", "m.2m", "m.2m", "m.2m"]
    @test dataset.hall_number == 419
    @test isempty(dataset.choice)
    @test dataset.equivalent_atoms == [0, 0, 2, 2, 2, 2]
    @test dataset.wyckoffs == ['a', 'a', 'f', 'f', 'f', 'f']
    @test dataset.std_positions == [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.3, 0.3, 0.0],
        [0.7, 0.7, 0.0],
        [0.2, 0.8, 0.5],
        [0.8, 0.2, 0.5],
    ]
    @test dataset.origin_shift == [0.0, 0.0, 0.0]
    @test dataset.international_symbol == "P4_2/mnm"
    @test get_international(rutile, 1e-5) == "P4_2/mnm"
    @test dataset.std_lattice == Lattice([
        4.0 0.0 0.0
        0.0 4.0 0.0
        0.0 0.0 3.0
    ])
    @test dataset.std_types == [1, 1, 2, 2, 2, 2]  # 14, 14, 8, 8, 8, 8
    @test dataset.pointgroup_symbol == "4/mmm"
    @test dataset.hall_symbol == "-P 4n 2n"
    @test dataset.transformation_matrix == [
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]
    @test dataset.std_mapping_to_primitive == [0, 1, 2, 3, 4, 5]
    @test dataset.crystallographic_orbits == [0, 0, 2, 2, 2, 2]
    @test dataset.primitive_lattice == Lattice([
        0.0 4.0 0.0
        0.0 0.0 4.0
        3.0 0.0 0.0
    ])
    @test dataset.std_rotation_matrix == [
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]
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
        [0 1 0; -1 0 0; 0 0 1],
        [0 -1 0; 1 0 0; 0 0 -1],
        [-1 0 0; 0 -1 0; 0 0 1],
        [1 0 0; 0 1 0; 0 0 -1],
        [0 -1 0; 1 0 0; 0 0 1],
        [0 1 0; -1 0 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 -1],
        [-1 0 0; 0 1 0; 0 0 1],
        [0 -1 0; -1 0 0; 0 0 -1],
        [0 1 0; 1 0 0; 0 0 1],
        [-1 0 0; 0 1 0; 0 0 -1],
        [1 0 0; 0 -1 0; 0 0 1],
        [0 1 0; 1 0 0; 0 0 -1],
        [0 -1 0; -1 0 0; 0 0 1],
    ]
    @test all(
        map(zip(dataset.rotations, python_rotations)) do (rotation, python_rotation)
            rotation == python_rotation
        end,
    )
    @test get_symmetry(rutile) == (dataset.rotations, dataset.translations)
end

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
    types = [14, 14, 8, 8, 8, 8]
    distorted_rutile = Cell(lattice, positions, types)
    @testset "Test `get_dataset`" begin
        # These results are compared with Python's spglib results.
        dataset = get_dataset(distorted_rutile, 1e-5)
        @test dataset.mapping_to_primitive == [0, 1, 2, 3, 4, 5]
        @test dataset.international_symbol == "P1"
        @test get_international(distorted_rutile, 1e-5) == "P1"
        @test dataset.spacegroup_number == 1
        @test dataset.site_symmetry_symbols == ["1", "1", "1", "1", "1", "1"]
        @test dataset.hall_number == 1
        @test isempty(dataset.choice)
        @test dataset.equivalent_atoms == [0, 1, 2, 3, 4, 5]
        @test only(dataset.rotations) == [
            1 0 0
            0 1 0
            0 0 1
        ]
        @test dataset.wyckoffs == ['a', 'a', 'a', 'a', 'a', 'a']
        @test dataset.primitive_lattice ==
            Lattice([[0.0, 0.0, 3.0], [3.97, 0.0, 0.0], [0.0, 4.03, 0.0]])
        @test dataset.std_positions == [
            [0.0, 0.0, 0.0],
            [0.5, 0.5001, 0.5],
            [0.0, 0.3, 0.3],
            [0.002, 0.7, 0.7],
            [0.5, 0.2, 0.8],
            [0.5, 0.8, 0.2],
        ]
        @test dataset.origin_shift == [0.0, 0.0, 0.0]
        @test dataset.std_lattice ≈ Lattice([
            [3, 0, 0], [2.4309239e-16, 3.97, 0], [2.4676633e-16, 2.4676633e-16, 4.03]
        ])
        @test dataset.std_types == [1, 1, 2, 2, 2, 2]
        @test dataset.std_rotation_matrix == [
            0.0 1.0 0.0
            0.0 0.0 1.0
            1.0 0.0 0.0
        ]
        @test dataset.translations == [[0.0, 0.0, 0.0]]
        @test dataset.pointgroup_symbol == "1"
        @test dataset.hall_symbol == "P 1"
        @test dataset.transformation_matrix == [
            0.0 1.0 0.0
            0.0 0.0 1.0
            1.0 0.0 0.0
        ]
        @test dataset.std_mapping_to_primitive == [0, 1, 2, 3, 4, 5]
        @test dataset.crystallographic_orbits == [0, 1, 2, 3, 4, 5]
        @test get_symmetry(distorted_rutile) == (dataset.rotations, dataset.translations)
    end
end

# From https://github.com/unkcpz/LibSymspg.jl/blob/53d2f6d/test/test_api.jl#L14-L32
@testset "Test `get_dataset` for a P-3m1 crystal" begin
    lattice = [
        4.0 2.0 0.0
        0.0 3.4641 0.0
        0.0 0.0 12.0
    ]
    positions = [[0.0, 0.0, 0.0], [1 / 3, 1 / 3, 1 / 3]]
    types = [1, 1]
    cell = Cell(lattice, positions, types)
    dataset = get_dataset(cell, 1e-3)
    # Compared results with Python spglib
    @test dataset.spacegroup_number == 164
    @test dataset.hall_number == 456
    @test dataset.international_symbol == "P-3m1"
    @test get_international(cell, 1e-3) == "P-3m1"
    @test isempty(dataset.choice)
    @test dataset.transformation_matrix ≈ [
        1.0 1.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]
    @test dataset.origin_shift ≈ [1 / 3, 2 / 3, 1 / 3]
    @test dataset.n_operations == 12
    python_rotations = [
        [1 0 0; 0 1 0; 0 0 1],
        [-1 0 0; 0 -1 0; 0 0 -1],
        [-1 1 0; -1 0 0; 0 0 1],
        [1 -1 0; 1 0 0; 0 0 -1],
        [0 -1 0; 1 -1 0; 0 0 1],
        [0 1 0; -1 1 0; 0 0 -1],
        [1 0 0; 1 -1 0; 0 0 -1],
        [-1 0 0; -1 1 0; 0 0 1],
        [0 -1 0; -1 0 0; 0 0 -1],
        [0 1 0; 1 0 0; 0 0 1],
        [-1 1 0; 0 1 0; 0 0 -1],
        [1 -1 0; 0 -1 0; 0 0 1],
    ]
    @test all(
        map(zip(dataset.rotations, python_rotations)) do (rotation, python_rotation)
            rotation == python_rotation
        end,
    )
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
    @test size(dataset.rotations) == (12,)
    @test size(dataset.translations) == (12,)
    @test get_symmetry(cell) == (dataset.rotations, dataset.translations)
    @test dataset.wyckoffs == ['d', 'd']
    @test dataset.site_symmetry_symbols == ["3m.", "3m."]
    @test dataset.equivalent_atoms == [0, 0]
    @test dataset.crystallographic_orbits == [0, 0]
    @test dataset.primitive_lattice ==
        Lattice([[0.0, 3.4641, 0.0], [-4.0, 1.4641, 0.0], [0.0, 0.0, 12.0]])
    @test dataset.mapping_to_primitive == [0, 1]
    @test dataset.std_lattice ≈
        Lattice([[3.4641, 0, 0], [-1.4641, 4, 0], [7.34788079e-16, 1.05141362e-15, 12]]) # ??
    @test dataset.std_types == [1, 1]
    @test dataset.std_positions ≈
        [[0.66666667, 0.83333333, 0.66666667], [0.33333333, 0.16666667, 0.33333333]]
    @test dataset.std_rotation_matrix ≈ [
        0.50000017 0.8660253 0.0
        -0.8660253 0.50000017 0.0
        0.0 0.0 1.0
    ]  # ??
    @test dataset.std_mapping_to_primitive == [0, 1]
    @test dataset.pointgroup_symbol == "-3m"
end

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
    types = [14, 14, 14, 14, 14, 14, 14, 14]
    silicon = Cell(lattice, positions, types)
    dataset = get_dataset(silicon, 1e-5)
    # Compared with Python results
    @test dataset.spacegroup_number == 227
    @test dataset.hall_number == 525
    @test dataset.international_symbol == "Fd-3m"
    @test dataset.hall_symbol == "F 4d 2 3 -1d"
    @test get_international(silicon, 1e-5) == "Fd-3m"
    @test dataset.choice == "1"
    @test dataset.transformation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
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
    ]  # Compared with Python `spglib`
    @test get_symmetry(silicon) == (dataset.rotations, dataset.translations)
    @test dataset.wyckoffs == ['a', 'a', 'a', 'a', 'a', 'a', 'a', 'a']
    @test dataset.site_symmetry_symbols ==
        ["-43m", "-43m", "-43m", "-43m", "-43m", "-43m", "-43m", "-43m"]
    @test dataset.equivalent_atoms == [0, 0, 0, 0, 0, 0, 0, 0]
    @test dataset.crystallographic_orbits == [0, 0, 0, 0, 0, 0, 0, 0]
    @test dataset.primitive_lattice ==
        Lattice([[2.0, -2.0, 0.0], [-2.0, -0.0, -2.0], [2.0, 2.0, 0.0]])
    @test dataset.mapping_to_primitive == [0, 0, 0, 0, 1, 1, 1, 1]
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
    @test dataset.std_mapping_to_primitive == [0, 1, 0, 1, 0, 1, 0, 1]
    @test dataset.pointgroup_symbol == "m-3m"
    @testset "Test primitive silicon structure" begin
        lattice = ([[0, 2, 2], [2, 0, 2], [2, 2, 0]])
        positions = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]
        types = [14, 14]
        silicon_prim = Cell(lattice, positions, types)
        dataset_prim = get_dataset(silicon_prim, 1e-5)
        for f in (
            :spacegroup_number,
            :hall_number,
            :hall_symbol,
            :international_symbol,
            :choice,
            :std_lattice,
            :std_rotation_matrix,
            :std_mapping_to_primitive,
            :pointgroup_symbol,
        )
            @test getfield(dataset_prim, f) == getfield(dataset, f)
        end
        @test dataset_prim.transformation_matrix == [
            0.0 0.5 0.5
            0.5 0.0 0.5
            0.5 0.5 0.0
        ]
        @test dataset_prim.origin_shift == [0, 0, 1 / 2]
        @test size(dataset_prim.rotations) == (48,)
        @test dataset_prim.std_types == [14, 14, 14, 14, 14, 14, 14, 14] / 14
        @test dataset_prim.std_positions == [
            [0.0, 0.0, 0.5],
            [0.25, 0.75, 0.25],
            [0.0, 0.5, 0.0],
            [0.25, 0.25, 0.75],
            [0.5, 0.0, 0.0],
            [0.75, 0.75, 0.75],
            [0.5, 0.5, 0.5],
            [0.75, 0.25, 0.25],
        ]
        @test dataset_prim.translations == [
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
        @test get_symmetry(silicon_prim) ==
            (dataset_prim.rotations, dataset_prim.translations)
        @test dataset_prim.wyckoffs == ['b', 'b']
        @test dataset_prim.site_symmetry_symbols == ["-43m", "-43m"]
        @test dataset_prim.equivalent_atoms == [0, 0]
        @test dataset_prim.crystallographic_orbits == [0, 0]
        @test dataset_prim.primitive_lattice ==
            Lattice([[2.0, -2.0, 0.0], [-2.0, -0.0, -2.0], [2.0, 2.0, 0.0]])
        @test dataset_prim.mapping_to_primitive == [0, 1]
    end
end
