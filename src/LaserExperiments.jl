module LaserExperiments

hello_world() = "Hello, World!"

# laser beams
export LaserBeam, laserbeam
export wavelength, classical_nonlinearity_parameter, unit_wavevector, intensity, frequency, photon_energy, four_momentum

using Unitful
using LinearAlgebra
using StaticArrays
using QEDcore

abstract type AbstractBeam end

include("units.jl")

include("beams/laser.jl")


end
