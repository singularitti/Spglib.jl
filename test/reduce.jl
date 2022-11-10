# From https://github.com/unkcpz/LibSymspg.jl/blob/f342e72/test/runtests.jl#L83-L85
@testset "Niggli reduce" begin
    lattice = [
        4.0 0.0 0.0
        20.0 2.0 0.0
        0.0 0.0 12.0
    ]
    @test niggli_reduce(lattice, 1e-3) ≈ [
        0.0 4.0 0.0
        -2.0 0.0 0.0
        0.0 0.0 12.0
    ]
    cell = Cell(lattice, [[0., 0., 0.], [0.05, 0.05, 0.05]], [1, 1])
    rcell = niggli_reduce(cell)
    c1 = Ref(cell.lattice) .* cell.positions 
    c2 = Ref(rcell.lattice) .* rcell.positions 
    # Cartesian coordinates should remain the same
    @test c1 == c2
end
# From https://github.com/unkcpz/LibSymspg.jl/blob/f342e72/test/runtests.jl#L87-89
@testset "Delaunay reduce" begin
    lattice = [
        4.0 0.0 0.0
        20.0 2.0 0.0
        0.0 0.0 12.0
    ]
    @test delaunay_reduce(lattice, 1e-3) ≈ [
        0.0 -4.0 0.0
        2.0 0.0 0.0
        0.0 0.0 12.0
    ]
    cell = Cell(lattice, [[0., 0., 0.], [0.05, 0.05, 0.05]], [1, 1])
    rcell = niggli_reduce(cell)
    c1 = Ref(cell.lattice) .* cell.positions 
    c2 = Ref(rcell.lattice) .* rcell.positions 
    # Cartesian coordinates should remain the same
    @test c1 == c2
end
