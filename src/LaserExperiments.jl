module LaserExperiments

# constants
export ME_eV

# laser beams
export LaserBeam, laserbeam
export classical_nonlinearity_parameter, photon_energy, polar_angle
export wavelength, unit_wavevector, intensity, frequency, four_momentum

# electron beam
export ElectronBeam, electronbeam
export lorentz_factor, beta, momentum_magnitude, unit_direction, four_momentum, total_energy

using Unitful
using LinearAlgebra
using StaticArrays
using QEDcore

abstract type AbstractBeam end

include("units.jl")
include("constants.jl")

include("beams/utils.jl")
include("beams/laser.jl")
include("beams/electron.jl")


end
