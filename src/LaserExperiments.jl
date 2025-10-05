module LaserExperiments

# this line tests the test coverage

hello_world() = "Hello, World!"

# laser beams
export LaserBeam, laserbeam
export classical_nonlinearity_parameter, photon_energy, polar_angle
export wavelength, unit_wavevector, intensity, frequency, four_momentum

using Unitful
using LinearAlgebra
using StaticArrays
using QEDcore

abstract type AbstractBeam end

include("units.jl")

include("beams/laser.jl")


end
