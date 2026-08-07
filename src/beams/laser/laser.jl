struct LaserSystem{T <: Real} <: AbstractLaser
    id::String
    base::GenericLaser{T}
    metadata::LaserSystemMetadata
end
