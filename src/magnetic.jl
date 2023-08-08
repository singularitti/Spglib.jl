export MagneticCell, MagneticDataset

"""
    Cell(lattice, positions, types, magmoms=zeros(length(types)))

The basic input data type of `Spglib`.

Lattice parameters `lattice` are given by a ``3×3`` matrix with floating point values,
where ``𝐚``, ``𝐛``, and ``𝐜`` are given as columns.
Fractional atomic positions `positions` are given
by a vector of ``N`` vectors with floating point values, where ``N`` is the number of atoms.
Numbers to distinguish atomic species `types` are given by a list of ``N`` integers.
The collinear polarizations `magmoms` only work with `get_symmetry` and are given
as a list of ``N`` floating point values, or a vector of vectors.
"""
@struct_hash_equal struct MagneticCell{L,P,T,M} <: AbstractCell
    lattice::Lattice{L}
    positions::Vector{MVector{3,P}}
    atoms::Vector{T}
    magmoms::M
end
function MagneticCell(lattice, positions, atoms, magmoms)
    if !(lattice isa Lattice)
        lattice = Lattice(lattice)
    end
    N = length(atoms)
    if positions isa AbstractMatrix
        P = eltype(positions)
        if size(positions) == (3, 3)
            error("ambiguous `positions` size 3×3! Use a vector of `Vector`s instead!")
        elseif size(positions) == (3, N)
            positions = collect(eachcol(positions))
        elseif size(positions) == (N, 3)
            positions = collect(eachrow(positions))
        else
            throw(
                DimensionMismatch(
                    "the `positions` has a different number of atoms from the `types`!"
                ),
            )
        end
    else  # positions isa AbstractVector or a Tuple
        P = eltype(Base.promote_typeof(positions...))
        positions = collect(map(MVector{3,P}, positions))
    end
    L, T, M = eltype(lattice), eltype(atoms), typeof(magmoms)
    return MagneticCell{L,P,T,M}(lattice, positions, atoms, magmoms)
end
MagneticCell(cell::Cell, magmoms) =
    MagneticCell(cell.lattice, cell.positions, cell.atoms, magmoms)

natoms(cell::MagneticCell) = length(cell.atoms)

atomtypes(cell::MagneticCell) = unique(cell.atoms)

"""
    Lattice(cell::MagneticCell)

Get the lattice of a `MagneticCell`.
"""
Lattice(cell::MagneticCell) = cell.lattice

function _expand_cell(cell::MagneticCell)
    lattice, positions, types, magmoms = cell.lattice,
    cell.positions, cell.atoms,
    cell.magmoms
    # Reference: https://github.com/mdavezac/spglib.jl/blob/master/src/spglib.jl#L32-L35 and https://github.com/spglib/spglib/blob/444e061/python/spglib/spglib.py#L953-L975
    clattice = Base.cconvert(Matrix{Cdouble}, transpose(lattice))
    cpositions = Base.cconvert(Matrix{Cdouble}, reduce(hcat, positions))
    ctypes = Cint[findfirst(isequal(u), unique(types)) for u in types]
    magmoms = Base.cconvert(Vector{Cdouble}, magmoms)
    return clattice, cpositions, ctypes, magmoms
end

function get_symmetry_with_collinear_spin(cell::MagneticCell, symprec=1e-5)
    lattice, positions, atoms, magmoms = _expand_cell(cell)
    n = length(cell.magmoms)
    # See https://github.com/spglib/spglib/blob/42527b0/python/spglib/spglib.py#L270
    max_size = 96n  # 96 = 48 × 2 since we have spins
    rotations = Array{Cint,3}(undef, 3, 3, max_size)
    translations = Matrix{Cdouble}(undef, 3, max_size)
    equivalent_atoms = Vector{Cint}(undef, n)
    nsym = @ccall libspglib.spg_get_symmetry_with_collinear_spin(
        rotations::Ptr{Cint},
        translations::Ptr{Cdouble},
        equivalent_atoms::Ptr{Cint},
        max_size::Cint,
        lattice::Ptr{Cdouble},
        positions::Ptr{Cdouble},
        atoms::Ptr{Cint},
        spins::Ptr{Cdouble},
        n::Cint,
        symprec::Cdouble,
    )::Cint
    check_error()
    rotations, translations = map(
        SMatrix{3,3,Int32,9}, eachslice(rotations[:, :, 1:nsym]; dims=3)
    ),
    map(SVector{3,Float64}, eachcol(translations[:, 1:nsym]))
    return rotations[:, :, 1:nsym], translations[:, 1:nsym]
end
const get_magnetic_symmetry = spg_get_symmetry_with_collinear_spin

struct SpglibMagneticDataset
    uni_number::Cint
    msg_type::Cint
    hall_number::Cint
    tensor_rank::Cint
    n_operations::Cint
    rotations::Ptr{Cint}
    translations::Ptr{Cdouble}
    time_reversals::Ptr{Cint}
    n_atoms::Cint
    equivalent_atoms::Ptr{Cint}
    transformation_matrix::NTuple{9,Cdouble}
    origin_shift::NTuple{3,Cdouble}
    n_std_atoms::Cint
    std_lattice::NTuple{9,Cdouble}
    std_types::Ptr{Cint}
    std_positions::Ptr{Cdouble}
    std_tensors::Ptr{Cdouble}
    std_rotation_matrix::NTuple{9,Cdouble}
    primitive_lattice::NTuple{9,Cdouble}
end

struct MagneticDataset
    uni_number::Int32
    msg_type::Int32
    hall_number::Int32
    tensor_rank::Int32
    n_operations::Int32
    rotations::Vector{SMatrix{3,3,Int32,9}}
    translations::Vector{SVector{3,Float64}}
    time_reversals::Vector{Int32}
    n_atoms::Int32
    equivalent_atoms::Vector{Int32}
    transformation_matrix::SMatrix{3,3,Float64,9}
    origin_shift::SVector{3,Float64}
    n_std_atoms::Int32
    std_lattice::Lattice{Float64}
    std_types::Vector{Int32}
    std_positions::Vector{SVector{3,Float64}}
    std_tensors::Vector{Float64}
    std_rotation_matrix::SMatrix{3,3,Float64,9}
    primitive_lattice::Lattice{Float64}
end

function Base.convert(::Type{MagneticDataset}, dataset::SpglibMagneticDataset)
    rotations = [
        _convert(SMatrix{3,3,Int32}, unsafe_load(dataset.rotations, i)) for
        i in Base.OneTo(dataset.n_operations)
    ]
    translations = [
        SVector{3}(unsafe_load(dataset.translations, i)) for
        i in Base.OneTo(dataset.n_operations)
    ]
    time_reversals = unsafe_wrap(
        Vector{Int32}, dataset.time_reversals, dataset.n_operations
    )
    equivalent_atoms = unsafe_wrap(Vector{Int32}, dataset.equivalent_atoms, dataset.n_atoms)
    std_lattice = Lattice(transpose(_convert(SMatrix{3,3,Float64}, dataset.std_lattice)))
    std_types = unsafe_wrap(Vector{Int32}, dataset.std_types, dataset.n_std_atoms)
    std_positions = [
        SVector{3}(unsafe_load(dataset.std_positions, i)) for
        i in Base.OneTo(dataset.n_std_atoms)
    ]
    std_tensors = unsafe_wrap(Vector{Float64}, dataset.std_tensors, dataset.n_std_atoms)
    return MagneticDataset(
        dataset.uni_number,
        dataset.msg_type,
        dataset.hall_number,
        dataset.tensor_rank,
        dataset.n_operations,
        rotations,
        translations,
        time_reversals,
        dataset.n_atoms,
        equivalent_atoms,
        dataset.transformation_matrix,
        dataset.origin_shift,
        dataset.n_std_atoms,
        std_lattice,
        std_types,
        std_positions,
        std_tensors,
        dataset.std_rotation_matrix,
        dataset.primitive_lattice,
    )
end
