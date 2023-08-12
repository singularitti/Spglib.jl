# See https://github.com/singularitti/Spglib.jl/issues/91#issuecomment-1206106977
@testset "Test example given by Jae-Mo Lihm (@jaemolihm)" begin
    lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    positions = [[-0.1, -0.1, -0.1], [0.1, 0.1, 0.1]]
    atoms = [1, 1]
    magmoms = [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]
    cell = Cell(lattice, positions, atoms, magmoms)
    rotations, translations = get_symmetry(cell, 1e-5)
    @test size(rotations) == size(translations) == (12,)
    rotations, translations, equivalent_atoms = get_magnetic_symmetry(cell, 1e-5)
    @test size(rotations) == size(translations) == (4,)
end
