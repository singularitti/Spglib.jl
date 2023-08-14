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
