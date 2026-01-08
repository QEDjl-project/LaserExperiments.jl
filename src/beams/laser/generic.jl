abstract type AbstractLaser end

struct GenericLaser{T <: Real} <: AbstractLaser
    name::String
    medium::LaserMedium
    classification::LaserClassification

    fundamental::FundamentalLaserParameters{T}
    beam::BeamParameters{T}
    optics::SystemOpticsParameters{T}
end
