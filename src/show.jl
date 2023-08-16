function Base.show(io::IO, ::MIME"text/plain", cell::SpglibCell)
    summary(io, cell)
    println(io)
    println(io, " lattice:")
    for row in eachrow(cell.lattice)
        println(io, "   ", join(row, "  "))
    end
    num_atom = natoms(cell)
    println(io, " $num_atom atomic positions:")
    for position in cell.positions
        println(io, "   ", join(position, "  "))
    end
    println(io, " $num_atom atoms:")
    println(io, "   ", join(cell.atoms, "  "))
    if !isempty(cell.magmoms)
        println(io, " $num_atom magmoms:")
        if eltype(cell.magmoms) <: AbstractArray
            for magmom in cell.magmoms
                println(io, "   ", magmom)
            end
        else
            println(io, "   ", join(cell.magmoms, "  "))
        end
    end
    return nothing
end
function Base.show(io::IO, ::MIME"text/plain", spgtype::SpacegroupType)
    summary(io, spgtype)
    println(io)
    for name in propertynames(spgtype)
        println(io, "  ", name, ": ", getfield(spgtype, name))
    end
    return nothing
end
function Base.show(io::IO, ::MIME"text/plain", dataset::Dataset)
    summary(io, dataset)
    println(io)
    for name in propertynames(dataset)
        println(io, "  ", name, ": ", getfield(dataset, name))
    end
    return nothing
end
