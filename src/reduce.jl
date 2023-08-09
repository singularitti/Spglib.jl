"""
    niggli_reduce(lattice::AbstractMatrix, symprec=1e-5)
    niggli_reduce(cell::Cell, symprec=1e-5)

Apply Niggli reduction to input basis vectors `lattice`.
"""
function niggli_reduce(lattice::Lattice, symprec=1e-5)
    clattice = transpose(convert(Matrix{Cdouble}, lattice))
    @ccall libsymspg.spg_niggli_reduce(clattice::Ptr{Cdouble}, symprec::Cdouble)::Cint
    check_error()
    return Lattice(transpose(clattice))
end
function niggli_reduce(cell::Cell, symprec=1e-5)
    lattice = Lattice(cell)
    clattice = niggli_reduce(lattice, symprec)
    # Keeping cartesian coordinates, see #106
    recip = inv(clattice) * cell.lattice
    new_frac_pos = [recip * pos for pos in cell.positions]
    return Cell(clattice, new_frac_pos, cell.atoms)
end

"""
    delaunay_reduce(lattice::AbstractMatrix, symprec=1e-5)
    delaunay_reduce(cell::Cell, symprec=1e-5)

Apply Delaunay reduction to input basis vectors `lattice`.
"""
function delaunay_reduce(lattice::Lattice, symprec=1e-5)
    clattice = transpose(convert(Matrix{Cdouble}, lattice))
    @ccall libsymspg.spg_delaunay_reduce(clattice::Ptr{Cdouble}, symprec::Cdouble)::Cint
    check_error()
    return Lattice(transpose(clattice))
end
function delaunay_reduce(cell::Cell, symprec=1e-5)
    lattice = Lattice(cell)
    clattice = delaunay_reduce(lattice, symprec)
    # Keeping cartesian coordinates, see #106
    recip = inv(clattice) * cell.lattice
    new_frac_pos = [recip * pos for pos in cell.positions]
    return Cell(clattice, new_frac_pos, cell.atoms)
end
