# TODO: use internal unit system to avoid these constants

abstract type AbstractLaserBeam <: AbstractBeam end

struct LaserBeam{T} <: AbstractLaserBeam
    a0::T
    wavelength::Quantity{T,Unitful.𝐋}
    unit_wavevector::SVector{3,T}
end


# Convenience constructor (keyword arguments)
function laserbeam(; a0=nothing, intensity=nothing, wavelength=nothing, frequency=nothing,
                   photon_energy=nothing, unit_wavevector= [0.0,0.0,1.0])

    # normalize direction
    k = normalize(SVector{3}(unit_wavevector))

    # infer missing quantities
    if wavelength === nothing && frequency !== nothing
        wavelength = u"c0" / frequency
    elseif wavelength === nothing && photon_energy!== nothing
        wavelength = (Unitful.ħ*Unitful.c0/photon_energy) |> u"µm"
    end

    if a0 === nothing && intensity !== nothing && wavelength !== nothing
        a0 = compute_a0(intensity, wavelength)
    end


    return LaserBeam(a0, wavelength, k)
end

# Derived quantities
# TODO: add docstrings with respective references
intensity(L::LaserBeam) = compute_intensity(L.a0, L.wavelength)
frequency(L::LaserBeam) = Unitful.c0 / L.wavelength
photon_energy(L::LaserBeam) = Unitful.ħ * frequency(L)
four_momentum(L::LaserBeam) = SVector(photon_energy(L)/Unitful.c0, photon_energy(L)/Unitful.c0 .* L.unit_wavevector)

function compute_a0(intensity::Quantity{T},wavelength::Quantity{T,Unitful.𝐋}) where T
    I = ustrip(uconvert(u"W"/u"cm^2",intensity)) / 10^22
    lambda = ustrip(uconvert(u"µm",wavelength))
    return 60*lambda*sqrt(I)
end

function compute_intensity(a0::T,wavelength::Quantity{T,Unitful.𝐋}) where T
    lambda = ustrip(uconvert(u"µm",wavelength))
    return (a0/(60*lambda))^2
end
