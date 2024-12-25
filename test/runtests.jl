using QuadPol
using Test

@testset "QuadPol.jl" begin
    include("test_calculations.jl")
    include("test_HR4796.jl")
end
