using Spglib
using Test

@testset "Spglib.jl" begin
    # Write your own tests here.
    include("symmetry.jl")
    include("standardize.jl")
    include("reciprocal.jl")
    include("reduce.jl")
    include("dataset.jl")
end
