module LaserExperiments

# constants
export ME_eV

# laser beams
export LaserParameters, laserparameter
export classical_nonlinearity_parameter, photon_energy, pulse_duration, pulse_energy, repetition_rate, coherence_length, polarization_ratio
#export wavelength, unit_wavevector, intensity, frequency, four_momentum

# electron beam
export ElectronBeam, electronbeam
export lorentz_factor, beta, momentum_magnitude, unit_direction, four_momentum, total_energy, polar_angle

using Unitful
using LinearAlgebra
using StaticArrays
using QEDcore

abstract type AbstractBeam end

include("units.jl")
include("constants.jl")

include("beams/utils.jl")
include("beams/laser/utils.jl")
include("beams/laser/parameter.jl")
include("beams/electron.jl")


end
