"""
    niggli_reduce(lattice::AbstractMatrix, symprec=1e-5)
    niggli_reduce(cell::Cell, symprec=1e-5)

Apply Niggli reduction to input basis vectors `lattice`.

The transformation from original basis vectors 
``\\begin{bmatrix} \\mathbf{a} & \\mathbf{b} & \\mathbf{c} \\end{bmatrix}``
to final basis vectors
``\\begin{bmatrix} \\mathbf{a}' & \\mathbf{b}' & \\mathbf{c}' \\end{bmatrix}``
is achieved by linear
combination of basis vectors with integer coefficients without
rotating coordinates. Therefore the transformation matrix is obtained
by

```math
\\mathbf{P} = \\begin{bmatrix} \\mathbf{a} & \\mathbf{b} & \\mathbf{c} \\end{bmatrix}
\\bigl(\\begin{bmatrix} \\mathbf{a}' & \\mathbf{b}' & \\mathbf{c}' \\end{bmatrix}\\bigr)^{-1}
```

and the matrix elements have to be almost integers.

See also [Computing rigid rotation introduced by idealization](@ref).
"""
function niggli_reduce(lattice::Lattice, symprec=1e-5)
    niggli_lattice = Base.cconvert(Matrix{Cdouble}, transpose(lattice))  # `transpose` must before `cconvert`!
    @ccall libsymspg.spg_niggli_reduce(niggli_lattice::Ptr{Cdouble}, symprec::Cdouble)::Cint
    check_error()
    return Lattice(transpose(niggli_lattice))
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

The transformation from original basis vectors 
``\\begin{bmatrix} \\mathbf{a} & \\mathbf{b} & \\mathbf{c} \\end{bmatrix}``
to final basis vectors
``\\begin{bmatrix} \\mathbf{a}' & \\mathbf{b}' & \\mathbf{c}' \\end{bmatrix}``
is achieved by linear
combination of basis vectors with integer coefficients without
rotating coordinates. Therefore the transformation matrix is obtained
by

```math
\\mathbf{P} = \\begin{bmatrix} \\mathbf{a} & \\mathbf{b} & \\mathbf{c} \\end{bmatrix}
\\bigl(\\begin{bmatrix} \\mathbf{a}' & \\mathbf{b}' & \\mathbf{c}' \\end{bmatrix}\\bigr)^{-1}
```

and the matrix elements have to be almost integers.

See also [Computing rigid rotation introduced by idealization](@ref).
"""
function delaunay_reduce(lattice::Lattice, symprec=1e-5)
    delaunay_lattice = Base.cconvert(Matrix{Cdouble}, transpose(lattice))  # `transpose` must before `cconvert`!
    @ccall libsymspg.spg_delaunay_reduce(
        delaunay_lattice::Ptr{Cdouble}, symprec::Cdouble
    )::Cint
    check_error()
    return Lattice(transpose(delaunay_lattice))
end
function delaunay_reduce(cell::Cell, symprec=1e-5)
    lattice = Lattice(cell)
    clattice = delaunay_reduce(lattice, symprec)
    # Keeping cartesian coordinates, see #106
    recip = inv(clattice) * cell.lattice
    new_frac_pos = [recip * pos for pos in cell.positions]
    return Cell(clattice, new_frac_pos, cell.atoms)
end
