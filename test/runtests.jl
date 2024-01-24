using Spglib
using Test

@testset "Spglib.jl" begin
    # Write your own tests here.
    include("misc.jl")
    include("magnetic.jl")
    include("standardize.jl")
    include("reciprocal.jl")
    include("reduce.jl")
    include("symmetry.jl")
    include("error.jl")
end
