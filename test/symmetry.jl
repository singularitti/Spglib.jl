function type2dict(dt)
    dict = Dict{Symbol,Any}()
    for n in propertynames(dt)
        dict[n] = getproperty(dt, n)
    end
    return dict
end

@testset "Test `get_spacegroup_type`" begin
    # Adapted from https://github.com/unkcpz/LibSymspg.jl/blob/53d2f6d/test/test_api.jl#L7-L12
    @test type2dict(get_spacegroup_type(101)) == Dict(
        :number => 15,
        :international_short => "C2/c",
        :international_full => "I 1 1 2/a",
        :international => "C 2/c = I 1 1 2/a",
        :schoenflies => "C2h^6",
        :hall_number => 101,
        :hall_symbol => "-I 2a",
        :choice => "-c3",
        :pointgroup_international => "2/m",
        :pointgroup_schoenflies => "C2h",
        :arithmetic_crystal_class_number => 8,
        :arithmetic_crystal_class_symbol => "2/mC",
    )
    @test type2dict(get_spacegroup_type(419)) == Dict(
        :number => 136,
        :international_short => "P4_2/mnm",
        :international_full => "P 4_2/m 2_1/n 2/m",
        :international => "P 4_2/m n m",
        :schoenflies => "D4h^14",
        :hall_number => 419,
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
        :hall_number => 1,
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
        :hall_number => 525,
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
        :hall_number => 485,
        :hall_symbol => "-P 6 2",
        :choice => "",
        :pointgroup_schoenflies => "D6h",
        :pointgroup_international => "6/mmm",
        :arithmetic_crystal_class_number => 58,
        :arithmetic_crystal_class_symbol => "6/mmm",
    )
end

# @testset "Test MgB2 structure" begin
#     a = 3.07
#     c = 3.52
#     lattice = [
#         a 0 0
#         -a/2 a/2*sqrt(3) 0
#         0 0 c
#     ]
#     positions = [
#         0 0 0
#         1.0/3 2.0/3 0.5
#         2.0/3 1.0/3 0.5
#     ]
#     types = [12, 5, 5]
#     MgB2 = Cell(lattice, positions, types)
# end

# From https://github.com/unkcpz/LibSymspg.jl/blob/53d2f6d/test/test_api.jl#L34-L77
@testset "Get symmetry operations" begin
    @testset "Normal symmetry" begin
        lattice = [
            4.0 0.0 0.0
            0.0 4.0 0.0
            0.0 0.0 4.0
        ]
        positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        types = [1, 1]
        cell = Cell(lattice, positions, types)
        num_atom = length(types)
        max_size = num_atom * 48
        rotation = Array{Cint,3}(undef, 3, 3, max_size)
        translation = Array{Float64,2}(undef, 3, max_size)
        get_symmetry!(rotation, translation, cell, 1e-5)
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
        positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        types = [1, 1]
        @testset "Test ferromagnetism" begin
            magmoms = [1.0, 1.0]
            cell = Cell(lattice, positions, types, magmoms)
            rotation, translation, equivalent_atoms = get_symmetry_with_collinear_spin(
                cell, 1e-5
            )
            @test size(rotation) == (3, 3, 96)
            @test equivalent_atoms == [0, 0]
        end
        @testset "Test antiferromagnetism" begin
            magmoms = [1.0, -1.0]
            cell = Cell(lattice, positions, types, magmoms)
            rotation, translation, equivalent_atoms = get_symmetry_with_collinear_spin(
                cell, 1e-5
            )
            @test size(rotation) == (3, 3, 96)
            @test equivalent_atoms == [0, 0]
        end
        @testset "Test broken magmoms" begin
            magmoms = [1.0, 2.0]
            cell = Cell(lattice, positions, types, magmoms)
            rotation, translation, equivalent_atoms = get_symmetry_with_collinear_spin(
                cell, 1e-5
            )
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
        positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        types = [1, 1]
        cell = Cell(lattice, positions, types)
        @test get_multiplicity(cell, 1e-5) == 96
    end
end
