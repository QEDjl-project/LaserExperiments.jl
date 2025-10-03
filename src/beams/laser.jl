# TODO: use internal unit system to avoid these constants

abstract type AbstractLaserBeam <: AbstractBeam end

struct LaserBeam{T} <: AbstractLaserBeam
    a0::T
    wavelength::Quantity{T, Unitful.𝐋}
    unit_wavevector::SVector{3, T}

    # automatically convert wavelength to µm
    function LaserBeam(a0::T, wavelength::Quantity{T, Unitful.𝐋}, unit_wavevector::SVector{3, T}) where {T <: Real}
        return new{T}(
            a0,
            uconvert(u"µm", wavelength),
            unit_wavevector
        )
    end

end


# Convenience constructor (keyword arguments)
function laserbeam(;
        a0 = nothing, intensity = nothing, wavelength = nothing, frequency = nothing,
        photon_energy = nothing, unit_wavevector = [0.0, 0.0, 1.0]
    )

    # normalize direction
    k = normalize(SVector{3}(unit_wavevector))

    # infer missing quantities
    if wavelength === nothing && frequency !== nothing
        wavelength = Unitful.c0 / frequency |> u"µm"
    elseif wavelength === nothing && photon_energy !== nothing
        wavelength = (Unitful.h * Unitful.c0 / photon_energy) |> u"µm"
    end

    if a0 === nothing && intensity !== nothing && wavelength !== nothing
        a0 = compute_a0(intensity, wavelength)
    end

    return LaserBeam(a0, wavelength, k)
end

# accessors
wavelength(L::LaserBeam) = L.wavelength
classical_nonlinearity_parameter(L::LaserBeam) = L.a0
unit_wavevector(L::LaserBeam) = L.unit_wavevector

# Derived quantities
# TODO:
# - add docstrings with respective references
# - check: frequency vs cicular frequency (what is used in Laser physics?)
intensity(L::LaserBeam) = compute_intensity(L.a0, L.wavelength)
frequency(L::LaserBeam) = uconvert(u"Hz", Unitful.c0 / L.wavelength)
photon_energy(L::LaserBeam) = mynatural(frequency(L))
function four_momentum(L::LaserBeam{T}) where {T}
    omega = photon_energy(L)
    unit_three_mom = unit_wavevector(L)
    return omega * SFourMomentum(one(T), unit_three_mom...)
end

function compute_a0(intensity::Quantity{T}, wavelength::Quantity{T, Unitful.𝐋}) where {T}
    I = ustrip(uconvert(u"W/cm^2", intensity)) / 1.0e22
    lambda = ustrip(uconvert(u"µm", wavelength))
    return 60 * lambda * sqrt(I)
end

function compute_intensity(a0::T, wavelength::Quantity{T, Unitful.𝐋}) where {T}
    lambda = ustrip(uconvert(u"µm", wavelength))
    return 1u"W/cm^2" * (a0 / (60 * lambda))^2 * 1.0e22
end
