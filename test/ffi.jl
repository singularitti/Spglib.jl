using Test

using Spglib

function type2dict(dt)
    di = Dict{Symbol,Any}()
    for n in propertynames(dt)
        di[n] = getproperty(dt, n)
    end
    di
end

@testset "Test `get_spacegroup_type`" begin
    # Adapted from https://github.com/unkcpz/LibSymspg.jl/blob/53d2f6d/test/test_api.jl#L7-L12
    spacegroup_type = get_spacegroup_type(101)
    @test spacegroup_type.number == 15
    @test spacegroup_type.hall_symbol == "-I 2a"
    @test spacegroup_type.arithmetic_crystal_class_symbol == "2/mC"
    # These results are compared with Python's spglib results.
    @test type2dict(get_spacegroup_type(419)) == Dict(
        :number => 136,
        :international_short => "P4_2/mnm",
        :international_full => "P 4_2/m 2_1/n 2/m",
        :international => "P 4_2/m n m",
        :schoenflies => "D4h^14",
        :hall_symbol => "-P 4n 2n",
        :choice => "",
        :pointgroup_schoenflies => "D4h",
        :pointgroup_international => "4/mmm",
        :arithmetic_crystal_class_number => 36,
        :arithmetic_crystal_class_symbol => "4/mmmP",
    )
    @test type2dict(get_spacegroup_type(1)) == Dict(
        :number => 1,
        :international_short => "P1",
        :international_full => "P 1",
        :international => "P 1",
        :schoenflies => "C1^1",
        :hall_symbol => "P 1",
        :choice => "",
        :pointgroup_schoenflies => "C1",
        :pointgroup_international => "1",
        :arithmetic_crystal_class_number => 1,
        :arithmetic_crystal_class_symbol => "1P",
    )
    @test type2dict(get_spacegroup_type(525)) == Dict(
        :number => 227,
        :international_short => "Fd-3m",
        :international_full => "F 4_1/d -3 2/m",
        :international => "F d -3 m",
        :schoenflies => "Oh^7",
        :hall_symbol => "F 4d 2 3 -1d",
        :choice => "1",
        :pointgroup_schoenflies => "Oh",
        :pointgroup_international => "m-3m",
        :arithmetic_crystal_class_number => 72,
        :arithmetic_crystal_class_symbol => "m-3mF",
    )
    @test type2dict(get_spacegroup_type(485)) == Dict(
        :number => 191,
        :international_short => "P6/mmm",
        :international_full => "P 6/m 2/m 2/m",
        :international => "P 6/m m m",
        :schoenflies => "D6h^1",
        :hall_symbol => "-P 6 2",
        :choice => "",
        :pointgroup_schoenflies => "D6h",
        :pointgroup_international => "6/mmm",
        :arithmetic_crystal_class_number => 58,
        :arithmetic_crystal_class_symbol => "6/mmm",
    )
end

@testset "Test silicon structure" begin
    lattice = [
        4 0 0
        0 4 0
        0 0 4
    ]
    positions = [
        0.0 0.0 0.5 0.5 0.25 0.25 0.75 0.75
        0.0 0.5 0.0 0.5 0.25 0.75 0.25 0.75
        0.0 0.5 0.5 0.0 0.25 0.75 0.75 0.25
    ]
    types = [14, 14, 14, 14, 14, 14, 14, 14]
    silicon = Cell(lattice, positions, types)
    dataset = get_dataset(silicon, 1e-5)
    # Compared with Python results
    @test dataset.spacegroup_number == 227
    @test dataset.hall_number == 525
    @test dataset.international_symbol == "Fd-3m"
    @test dataset.choice == "1"
    @test dataset.transformation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.origin_shift == [0, 0, 0]
    @test size(dataset.rotations) == (3, 3, 192)
    @test dataset.translations == [
        0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.75 0.5 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0
        0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5
        0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.75 0.0 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5
    ]
    @test dataset.wyckoffs == ["a", "a", "a", "a", "a", "a", "a", "a"]
    @test dataset.site_symmetry_symbols ==
          ["-43m", "-43m", "-43m", "-43m", "-43m", "-43m", "-43m", "-43m"]
    @test dataset.equivalent_atoms == [0, 0, 0, 0, 0, 0, 0, 0]
    if get_version() >= v"1.15"
        @test dataset.crystallographic_orbits == [0, 0, 0, 0, 0, 0, 0, 0]
        @test dataset.primitive_lattice == [
            2.0 -2.0 2.0
            -2.0 -0.0 2.0
            0.0 -2.0 0.0
        ]
    end
    @test dataset.mapping_to_primitive == [0, 0, 0, 0, 1, 1, 1, 1]
    @test dataset.std_lattice == 4 * [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.std_types == [14, 14, 14, 14, 14, 14, 14, 14] / 14
    @test dataset.std_positions == [
        0.0 0.25 0.0 0.25 0.5 0.75 0.5 0.75
        0.0 0.75 0.5 0.25 0.0 0.75 0.5 0.25
        0.0 0.75 0.5 0.25 0.5 0.25 0.0 0.75
    ]
    @test dataset.std_rotation_matrix == [
        1 0 0
        0 1 0
        0 0 1
    ]
    @test dataset.std_mapping_to_primitive == [0, 1, 0, 1, 0, 1, 0, 1]
    @test dataset.pointgroup_symbol == "m-3m"
    @testset "Test primitive silicon structure" begin
        lattice = [
            0 2 2
            2 0 2
            2 2 0
        ]
        positions = [
            0.0 0.25
            0.0 0.25
            0.0 0.25
        ]
        types = [14, 14]
        silicon_prim = Cell(lattice, positions, types)
        dataset_prim = get_dataset(silicon_prim)
        for f in (
            :spacegroup_number,
            :hall_number,
            :international_symbol,
            :choice,
            :std_lattice,
            :std_rotation_matrix,
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
        @test size(dataset_prim.rotations) == (3, 3, 48)
        # @test dataset_prim.std_types == [14, 14] / 14
        # @test dataset_prim.std_positions == [
        #     0.0 0.25
        #     0.0 0.75
        #     0.5 0.25
        # ]
        @test dataset_prim.translations == [
            0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0
            0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0
            0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0 0.25 0.0
        ]
        @test dataset_prim.wyckoffs == ["b", "b"]
        @test dataset_prim.site_symmetry_symbols == ["-43m", "-43m"]
        @test dataset_prim.equivalent_atoms == [0, 0]
        if get_version() >= v"1.15"
            @test dataset_prim.crystallographic_orbits == [0, 0]
            @test dataset_prim.primitive_lattice == [
                2.0 -2.0 2.0
                -2.0 0.0 2.0
                0.0 -2.0 0.0
            ]
        end
        @test dataset_prim.mapping_to_primitive == [0, 1]
    end
end

@testset "Test MgB2 structure" begin
    a = 3.07
    c = 3.52
    lattice = [
        a 0 0
        -a/2 a/2*sqrt(3) 0
        0 0 c
    ]
    positions = [
        0 0 0
        1.0/3 2.0/3 0.5
        2.0/3 1.0/3 0.5
    ]
    types = [12, 5, 5]
    MgB2 = Cell(lattice, positions, types)
end

# From https://github.com/unkcpz/LibSymspg.jl/blob/53d2f6d/test/test_api.jl#L34-L77
@testset "Get symmetry operations" begin
    @testset "Normal symmetry" begin
        lattice = [
            4.0 0.0 0.0
            0.0 4.0 0.0
            0.0 0.0 4.0
        ]
        positions = [
            0.0 0.5
            0.0 0.5
            0.0 0.5
        ]
        types = [1, 1]
        cell = Cell(lattice, positions, types)
        num_atom = length(types)
        max_size = num_atom * 48
        rotation = Array{Cint,3}(undef, 3, 3, max_size)
        translation = Array{Float64,2}(undef, 3, max_size)
        _ = get_symmetry!(rotation, translation, max_size, cell, 1e-5)
        @test size(rotation) == (3, 3, 96)
        @test size(translation) == (3, 96)
        @test get_hall_number_from_symmetry(rotation, translation, max_size, 1e-5) == 529
    end
    # See https://github.com/spglib/spglib/blob/deb6695/python/test/test_collinear_spin.py#L18-L37
    @testset "Get symmetry with collinear spins" begin
        lattice = [
            4.0 0.0 0.0
            0.0 4.0 0.0
            0.0 0.0 4.0
        ]
        positions = [
            0.0 0.5
            0.0 0.5
            0.0 0.5
        ]
        types = [1, 1]
        @testset "Test ferromagnetism" begin
            magmoms = [1.0, 1.0]
            cell = Cell(lattice, positions, types, magmoms)
            rotation, translation, equivalent_atoms =
                get_symmetry_with_collinear_spin(cell, 1e-5)
            @test size(rotation) == (3, 3, 96)
            @test equivalent_atoms == [0, 0]
        end
        @testset "Test antiferromagnetism" begin
            magmoms = [1.0, -1.0]
            cell = Cell(lattice, positions, types, magmoms)
            rotation, translation, equivalent_atoms =
                get_symmetry_with_collinear_spin(cell, 1e-5)
            @test size(rotation) == (3, 3, 96)
            @test equivalent_atoms == [0, 0]
        end
        @testset "Test broken magmoms" begin
            magmoms = [1.0, 2.0]
            cell = Cell(lattice, positions, types, magmoms)
            rotation, translation, equivalent_atoms =
                get_symmetry_with_collinear_spin(cell, 1e-5)
            @test size(rotation) == (3, 3, 48)
            @test size(translation) == (3, 48)
            @test equivalent_atoms == [0, 1]
        end
    end

    @testset "Get multiplicity" begin
        lattice = [
            4.0 0.0 0.0
            0.0 4.0 0.0
            0.0 0.0 4.0
        ]
        positions = [
            0.0 0.5
            0.0 0.5
            0.0 0.5
        ]
        types = [1, 1]
        cell = Cell(lattice, positions, types)
        @test get_multiplicity(cell, 1e-5) == 96
    end
end
