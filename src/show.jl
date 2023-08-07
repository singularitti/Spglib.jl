function Base.show(io::IO, ::MIME"text/plain", cell::MagneticCell)
    summary(io, cell)
    println(io)
    println(io, " lattice:")
    for row in eachrow(cell.lattice)
        println(io, "  ", join(row, "  "))
    end
    N = natoms(cell)
    println(io, " $N atomic positions:")
    for position in cell.positions
        println(io, "  ", position)
    end
    println(io, " $N atoms:")
    println(io, "  ", cell.types)
    return nothing
end
