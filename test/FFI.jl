using Test

using SpgLib.DataModel
using SpgLib.FFI

@testset "Test rutile structure" begin
    lattice = [
        4  0  0
        0  4  0
        0  0  3
    ]
    positions = [
        0.0  0.0  0.0
        0.5  0.5  0.5
        0.3  0.3  0.0
        0.7  0.7  0.0
        0.2  0.8  0.5
        0.8  0.2  0.5
    ]
    numbers = [14, 14, 8, 8, 8, 8]
    rutile = Cell(lattice, positions, numbers)
    # get_symmetry(rutile; symprec = 1e-5)
end # testset