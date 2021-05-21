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
    dataset = get_dataset(cell)
    # Compared with documented results
    @test dataset.spacegroup_number == 64
    @test dataset.international_symbol == "Cmce"
    @test dataset.hall_number == 304  # Compared with Python results
    @test dataset.hall_symbol == "-C 2bc 2"  # Compared with Python results
    @test dataset.transformation_matrix == [  # Compared with documented results
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.std_lattice ≈ [  # Compared with documented results
        7.17851431 0 0
        0 3.99943947 0
        0 0 8.57154746
    ]
    # Compared with Python results
    @test isempty(dataset.choice)
    @test dataset.origin_shift ≈ [5.55111512e-17, 0, 0]
    @test dataset.pointgroup_symbol == "mmm"
    @test dataset.std_types == [35, 35, 35, 35, 35, 35, 35, 35] / 35
    @test dataset.std_positions ≈ [
        0.0 0.0 0.0 0.0 0.5 0.5 0.5 0.5
        0.84688439 0.65311561 0.34688439 0.15311561 0.34688439 0.15311561 0.84688439 0.65311561
        0.1203133 0.6203133 0.3796867 0.8796867 0.1203133 0.6203133 0.3796867 0.8796867
    ]
    @test dataset.std_rotation_matrix ≈ [
        0.70710678 -0.70710678 0.0
        0.70710678 0.70710678 0.0
        0.0 0.0 1.0
    ]
    @test dataset.mapping_to_primitive == [0, 1, 2, 3, 0, 1, 2, 3]
    @test dataset.std_mapping_to_primitive == [0, 1, 2, 3, 0, 1, 2, 3]
    @test dataset.wyckoffs == ["f", "f", "f", "f", "f", "f", "f", "f"]
    @test dataset.site_symmetry_symbols ==
          ["m..", "m..", "m..", "m..", "m..", "m..", "m..", "m.."]
    @test dataset.equivalent_atoms == [0, 0, 0, 0, 0, 0, 0, 0]
    if get_version() >= v"1.15"
        @test dataset.crystallographic_orbits == [0, 0, 0, 0, 0, 0, 0, 0]
        @test dataset.primitive_lattice ≈ [
            2.82803077 1.12397269 0.0
            -2.82803077 3.95200346 0.0
            0.0 0.0 8.57154746
        ]
    end
    @test dataset.translations ≈ [
        0.0 -1.11022302e-16 -1.11022302e-16 0.0 0.0 -1.11022302e-16 -1.11022302e-16 0.0 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
        0.0 2.46519033e-32 0.5 0.5 0.0 2.46519033e-32 0.5 0.5 0.5 0.5 0.0 0.0 0.5 0.5 0.0 0.0
        0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.5
    ]
    rotations = dataset.rotations
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
    @test all(map(1:16) do i
        rotations[:, :, i] == python_rotations[i]
    end)
end
