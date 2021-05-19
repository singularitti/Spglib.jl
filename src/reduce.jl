"""
    niggli_reduce!(lattice::AbstractMatrix, symprec=1e-5)

Apply Niggli reduction to input basis vectors `lattice` and the reduced basis vectors are overwritten to `lattice`.
"""
function niggli_reduce!(lattice::AbstractMatrix, symprec = 1e-5)
    clattice = convert(Matrix{Cdouble}, lattice)
    exitcode = ccall(
        (:spg_niggli_reduce, libsymspg),
        Cint,
        (Ptr{Cdouble}, Cdouble),
        clattice,
        symprec,
    )
    iszero(exitcode) && error("Niggli reduce failed!")
    return clattice
end
function niggli_reduce!(cell::Cell, symprec = 1e-5)
    clattice = niggli_reduce!(cell.lattice, symprec)
    cell.lattice[:, :] = clattice
    return cell
end

"""
    delaunay_reduce!(lattice::AbstractMatrix, symprec=1e-5)

Apply Delaunay reduction to input basis vectors `lattice` and the reduced basis vectors are overwritten to `lattice`.
"""
function delaunay_reduce!(lattice::AbstractMatrix, symprec = 1e-5)
    clattice = convert(Matrix{Cdouble}, lattice)
    exitcode = ccall(
        (:spg_delaunay_reduce, libsymspg),
        Cint,
        (Ptr{Cdouble}, Cdouble),
        clattice,
        symprec,
    )
    iszero(exitcode) && error("Delaunay reduce failed!")
    return clattice
end
function delaunay_reduce!(cell::Cell, symprec = 1e-5)
    clattice = delaunay_reduce!(cell.lattice, symprec)
    cell.lattice[:, :] = clattice
    return cell
end
