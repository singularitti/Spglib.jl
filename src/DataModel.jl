"""
# module DataModel



# Examples

```jldoctest
julia>
```
"""
module DataModel

export Cell

struct Cell{L<:AbstractMatrix,P<:AbstractMatrix,N<:AbstractVector}
    lattice::L
    positions::P
    numbers::N
end

end
