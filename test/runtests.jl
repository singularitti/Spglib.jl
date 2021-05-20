using Spglib
using Test

@testset "Spglib.jl" begin
    # Write your own tests here.
    include("ffi.jl")
    include("standardize.jl")
end
