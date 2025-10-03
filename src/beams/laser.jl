# TODO: use internal unit system to avoid these constants

abstract type AbstractLaserBeam <: AbstractBeam end

"""
    struct LaserBeam{T} <: AbstractLaserBeam
    LaserBeam(a0,wavelength,unit_wavevector)

A parametric type representing a monochromatic plane-wave laser beam.
Stores the dimensionless field strength parameter `a0`, the wavelength (internally converted to µm),
and a normalized propagation direction (`unit_wavevector`).

# Fields
- `a0::T`: Dimensionless classical nonlinearity parameter
  (a measure of the laser’s normalized vector potential).
- `wavelength::Quantity{T,𝐋}`: Laser wavelength in µm.
- `unit_wavevector::SVector{3,T}`: Unit vector giving the propagation direction of the laser.

# Notes
The wavelength is automatically converted to µm for consistent scaling in auxiliary functions.
"""
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

"""
    laserbeam(; a0, intensity, wavelength, frequency, photon_energy, unit_wavevector=[0,0,1])

Convenience constructor for `LaserBeam`.

Automatically infers missing parameters where possible:
- If `wavelength` is not given, it is inferred from `frequency` or `photon_energy`.
- If `a0` is not given, it is computed from the provided `intensity` and `wavelength`.

# Keyword arguments
- `a0`: Classical nonlinearity parameter (dimensionless).
- `intensity`: Laser intensity as a physical quantity (e.g. `1e20u"W/cm^2"`).
- `wavelength`: Wavelength of the laser (any length unit).
- `frequency`: Frequency of the laser (Hz).
- `photon_energy`: Photon energy (e.g. eV, J).
- `unit_wavevector`: Direction vector (normalized internally).

# Returns
A fully constructed `LaserBeam` with consistent parameters.
"""
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

"""
    wavelength(L::LaserBeam)

Return the wavelength of the laser beam in ``\\mu\\mathrm{m}``.
"""
wavelength(L::LaserBeam) = L.wavelength

"""
    classical_nonlinearity_parameter(L::LaserBeam)

Return the dimensionless laser strength parameter ``a_0``.
"""
classical_nonlinearity_parameter(L::LaserBeam) = L.a0

"""
    unit_wavevector(L::LaserBeam)

Return the normalized propagation direction of the laser beam as a 3D vector.
"""
unit_wavevector(L::LaserBeam) = L.unit_wavevector

# Derived quantities
# TODO:
# - add docstrings with respective references
# - check: frequency vs cicular frequency (what is used in Laser physics?)

"""
    intensity(L::LaserBeam)

Compute the laser intensity from `a0` and `wavelength`.
Returned as a physical quantity in ``\\mathrm{W}/\\mathrm{cm}^2``.
"""
intensity(L::LaserBeam) = compute_intensity(L.a0, L.wavelength)

"""
    frequency(L::LaserBeam)

Return the frequency of the laser in Hz, defined as ``\\nu = c/\\lambda``.

Note: This uses the linear frequency (Hz), not the angular frequency (rad/s).
"""
frequency(L::LaserBeam) = uconvert(u"Hz", Unitful.c0 / L.wavelength)

"""
    photon_energy(L::LaserBeam)

Return the photon energy in eV associated with the laser frequency.
"""
photon_energy(L::LaserBeam) = mynatural(frequency(L))

"""
    four_momentum(L::LaserBeam)

Return the photon four-momentum associated with the laser beam:

``k^\\mu = (\\omega, \\vec{k})``

where ``\\omega`` is the photon energy in natural units and ``\\vec{k}`` points along the normalized propagation direction.
"""
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
