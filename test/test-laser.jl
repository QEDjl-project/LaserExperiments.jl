using Test
using SafeTestsets

begin
    @safetestset "parameter" begin
        include("laser/parameter.jl")
    end
end
