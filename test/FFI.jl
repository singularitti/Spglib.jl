using Test

function type2dict(dt)
    di = Dict{Symbol,Any}()
    for n in propertynames(dt)
        di[n] = getproperty(dt, n)
    end
    di
end

using Spglib.DataModel: Cell, Dataset, SpaceGroup
using Spglib.FFI

@testset "Test rutile structure" begin
    lattice = [
        4 0 0
        0 4 0
        0 0 3
    ]
    positions = [
        0.0 0.0 0.0
        0.5 0.5 0.5
        0.3 0.3 0.0
        0.7 0.7 0.0
        0.2 0.8 0.5
        0.8 0.2 0.5
    ]
    numbers = [14, 14, 8, 8, 8, 8]
    rutile = Cell(lattice, positions, numbers)
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
    # get_symmetry(rutile; symprec = 1e-5)
end # testset

@testset "Test distorted rutile structure" begin
    lattice = [
        3.97 0.00 0.00
        0.00 4.03 0.00
        0.00 0.00 3.00
    ]
    positions = [
        0 0 0
        0.5001 0.5 0.5
        0.3 0.3 0.0
        0.7 0.7 0.002
        0.2 0.8 0.5
        0.8 0.2 0.5
    ]
    numbers = [14, 14, 8, 8, 8, 8]
    distorted_rutile = Cell(lattice, positions, numbers)
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
end # testset

@testset "Test silicon structure" begin
    lattice = [
        4 0 0
        0 4 0
        0 0 4
    ]
    positions = [
        0 0 0
        0 0.5 0.5
        0.5 0 0.5
        0.5 0.5 0
        0.25 0.25 0.25
        0.25 0.75 0.75
        0.75 0.25 0.75
        0.75 0.75 0.25
    ]
    numbers = [14, 14, 14, 14, 14, 14, 14, 14]
    silicon = Cell(lattice, positions, numbers)
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
end

@testset "Test silicon_dist structure" begin
    lattice = [
        4.01 0 0
        0 4 0
        0 0 3.99
    ]
    positions = [
        0.001 0 0
        0 0.5 0.5
        0.5 0 0.5
        0.5 0.5 0
        0.25 0.25 0.251
        0.25 0.75 0.75
        0.75 0.25 0.75
        0.75 0.75 0.25
    ]
    numbers = [14, 14, 14, 14, 14, 14, 14, 14]
    silicon_dist = Cell(lattice, positions, numbers)
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
end

@testset "Test silicon_prim structure" begin
    lattice = [
        0 2 2
        2 0 2
        2 2 0
    ]
    positions = [
        0 0 0
        0.25 0.25 0.25
    ]
    numbers = [14, 14]
    silicon_prim = Cell(lattice, positions, numbers)
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
    numbers = [12, 5, 5]
    MgB2 = Cell(lattice, positions, numbers)
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
