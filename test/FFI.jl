using Test

using SpgLib
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
end # testset