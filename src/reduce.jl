"""
    niggli_reduce(lattice::AbstractMatrix, symprec=1e-5)
    niggli_reduce(cell::Cell, symprec=1e-5)

Apply Niggli reduction to input basis vectors `lattice`.
"""
function niggli_reduce(lattice::AbstractMatrix, symprec = 1e-5)
    clattice = convert(Matrix{Cdouble}, transpose(lattice))
    exitcode = ccall(
        (:spg_niggli_reduce, libsymspg),
        Cint,
        (Ptr{Cdouble}, Cdouble),
        clattice,
        symprec,
    )
    iszero(exitcode) && error("Niggli reduce failed!")
    return transpose(clattice)
end
function niggli_reduce(cell::Cell, symprec = 1e-5)
    lattice = cell.lattice
    clattice = niggli_reduce(lattice, symprec)
    return Cell(clattice, cell.positions, cell.types, cell.magmoms)
end

"""
    delaunay_reduce(lattice::AbstractMatrix, symprec=1e-5)
    delaunay_reduce(cell::Cell, symprec=1e-5)

Apply Delaunay reduction to input basis vectors `lattice`.
"""
function delaunay_reduce(lattice::AbstractMatrix, symprec = 1e-5)
    clattice = convert(Matrix{Cdouble}, transpose(lattice))
    exitcode = ccall(
        (:spg_delaunay_reduce, libsymspg),
        Cint,
        (Ptr{Cdouble}, Cdouble),
        clattice,
        symprec,
    )
    iszero(exitcode) && error("Delaunay reduce failed!")
    return transpose(clattice)
end
function delaunay_reduce(cell::Cell, symprec = 1e-5)
    lattice = cell.lattice
    clattice = delaunay_reduce(lattice, symprec)
    return Cell(clattice, cell.positions, cell.types, cell.magmoms)
end
