module LaserExperiments

hello_world() = "Hello, World!"

using Unitful
using StaticArrays

abstract type AbstractBeam end

include("beams/laser.jl")

end
