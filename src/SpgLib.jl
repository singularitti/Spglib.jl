module SpgLib

using Parameters: @with_kw

export Cell

@with_kw struct Cell{L<:AbstractMatrix,P<:AbstractMatrix,N<:AbstractVector,M<:Union{AbstractVector,Nothing}}
    lattice::L
    positions::P
    numbers::N
    magmoms::M = nothing
end
Cell(lattice, positions, numbers) = Cell(lattice, positions, numbers, nothing)

include("Wrapper/Wrapper.jl")
include("FFI.jl")

end # module
