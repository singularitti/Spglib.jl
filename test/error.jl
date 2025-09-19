using Spglib:
    SUCCESS,
    SPACEGROUP_SEARCH_FAILED,
    NIGGLI_FAILED,
    DELAUNAY_FAILED,
    CELL_STANDARDIZATION_FAILED,
    ATOMS_TOO_CLOSE,
    SpglibError,
    get_error_code

@testset "Test when no error occurs" begin
    lattice = [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]] * 5.4
    positions = [[0.875, 0.875, 0.875], [0.125, 0.125, 0.125]]
    atoms = [1, 1]
    cell = Cell(lattice, positions, atoms)
    get_dataset(cell, 1e-5)
    @test get_error_code() == SUCCESS
end

@testset "Test error types" begin
    lattice = Lattice([[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]] * 5.4)
    positions = [[0.875, 0.875, 0.875], [0.125, 0.125, 0.125]]
    atoms = [1, 1]
    cell = Cell(lattice, positions, atoms)
    @test_throws SpglibError get_dataset(cell, 1e-5)
    @test get_error_code() == SPACEGROUP_SEARCH_FAILED
    @test_throws SpglibError niggli_reduce(cell, 1e-5)
    @test get_error_code() == NIGGLI_FAILED
    @test_throws SpglibError delaunay_reduce(cell, 1e-5)
    @test get_error_code() == DELAUNAY_FAILED
    @test_throws SpglibError standardize_cell(cell, 1e-5)
    @test get_error_code() == CELL_STANDARDIZATION_FAILED
    @testset "Test when atoms are too close" begin
        positions = [[0.875, 0.875, 0.875], [0.875 + eps(), 0.875, 0.875]]
        cell = Cell(lattice, positions, atoms)
        @test_throws SpglibError get_dataset(cell, 1e-5)
        @test get_error_code() == ATOMS_TOO_CLOSE
    end
end
