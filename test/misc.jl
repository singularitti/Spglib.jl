@test get_version() == v"2.5.0"

@testset "Test constructors" begin
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
    @test Lattice(cell) == Lattice(
        [
            5.07597614744567 -2.8280307701821314 0.0
            5.07597614744567 2.8280307701821314 0.0
            0.0 0.0 8.57154746
        ],
    )
    @test_throws DimensionMismatch Cell(lattice, positions[1:3], atoms)
    @test_throws DimensionMismatch Cell(
        lattice, positions, atoms, fill(1.0, length(positions) + 1)
    )
    @test_throws DimensionMismatch Cell(
        lattice, positions, atoms, fill([-1, 1, 0], length(positions) + 1)
    )
    @test_throws DimensionMismatch Cell(lattice, positions, atoms, zeros(length(atoms), 3))
    @test_throws DimensionMismatch Cell(lattice, positions, atoms, zeros(length(atoms), 2))
    @test_throws MethodError Cell(lattice, positions, atoms, zeros(length(atoms), 1))
    @test_throws DimensionMismatch Cell(
        lattice, positions, atoms, fill([1], length(positions))
    )
end
