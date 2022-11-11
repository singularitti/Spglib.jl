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
    if exitcode == 0
        throw(SpglibError("Niggli reduce failed!"))
    end
    return transpose(clattice)
end
function niggli_reduce(cell::Cell, symprec = 1e-5)
    lattice = cell.lattice
    clattice = niggli_reduce(lattice, symprec)
    # Keeping cartesian coordinates, see #106
    recip = inv(clattice) * cell.lattice
    new_frac_pos = [recip * pos for pos in cell.positions]
    return Cell(clattice, new_frac_pos, cell.types, cell.magmoms)
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
    if exitcode == 0
        throw(SpglibError("Delaunay reduce failed!"))
    end
    return transpose(clattice)
end
function delaunay_reduce(cell::Cell, symprec = 1e-5)
    lattice = cell.lattice
    clattice = delaunay_reduce(lattice, symprec)
    # Keeping cartesian coordinates, see #106
    recip = inv(clattice) * cell.lattice
    new_frac_pos = [recip * pos for pos in cell.positions]
    return Cell(clattice, new_frac_pos, cell.types, cell.magmoms)
end
