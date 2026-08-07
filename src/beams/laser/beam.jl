abstract type AbstractBeamProfile end

struct GaussianBeam <: AbstractBeamProfile end
struct TopHatBeam <: AbstractBeamProfile end

struct BeamParameters{T <: Real}
    diameter_mm::T                 # 1/e^2 definition
    power_density_W_cm2::T
    energy_density_uJ_cm2::T
    profile::AbstractBeamProfile
    divergence_mrad::T
end

struct SystemOpticsParameters{T <: Real}
    spot_size_um::T
    working_distance_um::T
end
