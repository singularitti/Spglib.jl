"""
# module DataModel



# Examples

```jldoctest
julia>
```
"""
module DataModel

using Parameters

export Cell

@with_kw struct Cell{L<:AbstractMatrix,P<:AbstractMatrix,N<:AbstractVector,M<:Union{AbstractVector,Nothing}}
    lattice::L
    positions::P
    numbers::N
    magmoms::M = nothing
end
Cell(lattice, positions, numbers) = Cell(lattice, positions, numbers, nothing)

end
