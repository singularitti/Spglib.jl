# From https://github.com/spglib/spglib/blob/d8c39f6/example/python_api/example_full.py#L119-L122
@testset "Test Niggli reduction from Python example" begin
    𝐚 = [3, 0, 0]
    𝐛 = [-3.66666667, 3.68178701, 0]
    𝐜 = [-0.66666667, -1.3429469, 1.32364995]
    lattice = Lattice(𝐚, 𝐛, 𝐜)
    @test niggli_reduce(lattice, 1e-5) ≈ Lattice([
        [-0.66666667, -1.3429469, 1.32364995],
        [2.33333333, -1.3429469, 1.32364995],
        [0.99999999, 0.99589321, 2.6472999],
    ])
end

# From https://github.com/unkcpz/LibSymspg.jl/blob/f342e72/test/runtests.jl#L83-L85
@testset "Test Niggli reduction" begin
    lattice = Lattice([
        4 0 0
        20 2 0
        0 0 12
    ])  # Note in `LibSymspg.jl`, the lattice is transposed
    @test niggli_reduce(lattice, 1e-3) == Lattice([[0, -2, 0], [4, 0, 0], [0, 0, 12]])  # Compared also with Python results
    cell = Cell(lattice, [[0.0, 0.0, 0.0], [0.05, 0.05, 0.05]], [1, 1])
    rcell = niggli_reduce(cell)
    c1 = Ref(cell.lattice) .* cell.positions
    c2 = Ref(rcell.lattice) .* rcell.positions
    # Cartesian coordinates should remain the same
    @test c1 == c2
end

# From https://github.com/unkcpz/LibSymspg.jl/blob/f342e72/test/runtests.jl#L87-89
@testset "Test Delaunay reduction" begin
    lattice = Lattice([
        4 0 0
        20 2 0
        0 0 12
    ])  # Note in `LibSymspg.jl`, the lattice is transposed
    @test delaunay_reduce(lattice, 1e-3) == Lattice([
        0 -4 0
        2 0 0
        0 0 12
    ])
    cell = Cell(lattice, [[0.0, 0.0, 0.0], [0.05, 0.05, 0.05]], [1, 1])
    rcell = delaunay_reduce(cell)
    c1 = Ref(cell.lattice) .* cell.positions
    c2 = Ref(rcell.lattice) .* rcell.positions
    # Cartesian coordinates should remain the same
    @test c1 == c2
end
