@testset "Test examples from Python" begin
    # From https://github.com/spglib/spglib/blob/d8c39f6/example/python_api/example_full.py#L119-L122
    @testset "Unknown lattice" begin
        𝐚 = [3, 0, 0]
        𝐛 = [-3.66666667, 3.68178701, 0]
        𝐜 = [-0.66666667, -1.3429469, 1.32364995]
        lattice = Lattice(𝐚, 𝐛, 𝐜)
        @test niggli_reduce(lattice, 1e-5) ≈ Lattice([
            [-0.66666667, -1.3429469, 1.32364995],
            [2.33333333, -1.3429469, 1.32364995],
            [0.99999999, 0.99589321, 2.6472999],
        ])  # Compared with Python result
        @test delaunay_reduce(lattice, 1e-5) ≈ Lattice([
            [-0.66666667, -1.3429469, 1.32364995],
            [-2.33333333, 1.3429469, -1.32364995],
            [-0.99999999, -0.99589321, -2.6472999],
        ])  # Compared with Python result
    end
    # From https://spglib.github.io/spglib/definition.html#computing-rigid-rotation-introduced-by-idealization
    @testset "Test another lattice" begin
        lattice = Lattice([
            [5.0759761474456697, 5.0759761474456697, 0],
            [-2.8280307701821314, 2.8280307701821314, 0],
            [0, 0, 8.57154746],
        ])
        @test niggli_reduce(lattice, 1e-5) ≈ Lattice([
            [2.82803077, -2.82803077, 0.0],
            [-5.07597615, -5.07597615, 0.0],
            [0.0, 0.0, -8.57154746],
        ])  # Compared with Python result
        @test delaunay_reduce(lattice, 1e-5) ≈ Lattice([
            [2.82803077, -2.82803077, 0.0],
            [-5.07597615, -5.07597615, 0.0],
            [0.0, 0.0, -8.57154746],
        ])  # Compared with Python result
    end
end

@testset "Test examples from `LibSymspg.jl`" begin
    lattice = Lattice([
        4 0 0
        20 2 0
        0 0 12
    ])  # Note in `LibSymspg.jl`, the lattice is transposed
    # From https://github.com/unkcpz/LibSymspg.jl/blob/f342e72/test/runtests.jl#L83-L85
    @testset "Test Niggli reduction" begin
        @test niggli_reduce(lattice, 1e-3) == Lattice([[0, -2, 0], [4, 0, 0], [0, 0, 12]])  # Compared also with Python results
        cell = Cell(lattice, [[0.0, 0.0, 0.0], [0.05, 0.05, 0.05]], [1, 1])
        rcell = niggli_reduce(cell)
        c1 = Ref(cell.lattice) .* cell.positions
        c2 = Ref(rcell.lattice) .* rcell.positions
        # Cartesian coordinates should remain the same before and after basis vectors changes
        @test c1 == c2
    end
    # From https://github.com/unkcpz/LibSymspg.jl/blob/f342e72/test/runtests.jl#L87-89
    @testset "Test Delaunay reduction" begin
        @test delaunay_reduce(lattice, 1e-3) == Lattice([
            0 -4 0
            2 0 0
            0 0 12
        ])
        cell = Cell(lattice, [[0.0, 0.0, 0.0], [0.05, 0.05, 0.05]], [1, 1])
        rcell = delaunay_reduce(cell)
        c1 = Ref(cell.lattice) .* cell.positions
        c2 = Ref(rcell.lattice) .* rcell.positions
        # Cartesian coordinates should remain the same before and after basis vectors changes
        @test c1 == c2
    end
end

# From https://github.com/spglib/spglib/blob/4aa0806/test/functional/c/test_delaunay.cpp
@testset "Test an example from C" begin
    lattice = Lattice([
        1 0 0
        0 37 10
        0 11 3
    ])
    @test niggli_reduce(lattice, 1e-5) == Lattice([[0, 1, 0], [1, 0, 0], [0, 0, -1]])  # Compared with Python result
    @test delaunay_reduce(lattice, 1e-5) ≈ Lattice([[1, 0, 0], [0, -1, 0], [0, 0, -1]])  # Compared with Python result
end
