abstract type AbstractLaserBeam <: AbstractBeam end

"""

    struct LaserBeam{T} <: AbstractLaserBeam
    LaserBeam(a0, omega, alpha)

Representation of a monochromatic plane-wave laser beam.

## Fields
- `a0::T`: Dimensionless classical nonlinearity parameter (laser strength parameter).
- `omega::Quantity{T}`: Photon energy (in eV).
- `alpha::T`: Polar angle (in radians) defining the laser propagation direction.

## Construction

A new LaserBeam can be created via the outer constructor laserbeam.
The inner constructor performs sanity checks and unit normalization.

## Example

```
L = LaserBeam(1.0, 1.55u"eV", π/4)
```
"""
struct LaserBeam{T} <: AbstractLaserBeam
    a0::T
    omega::Quantity{T}
    alpha::T

    # automatically convert wavelength to µm
    function LaserBeam(a0::T, omega::Quantity{T}, alpha::T) where {T <: Real}
        zero(T) <= alpha <= T(pi) || throw(
            ArgumentError(
                "polar angle `alpha` must be between 0 and pi"
            )
        )
        a0 >= zero(T) || throw(
            ArgumentError(
                "classical nonlinearity parameter `a0` must be positive"
            )
        )
        ustrip(omega) > zero(T) || throw(
            ArgumentError(
                "value of photon energy `omega` must be positive"
            )
        )

        return new{T}(a0, mynatural(omega), alpha)
    end

end

"""
    laserbeam(; a0, intensity, wavelength, frequency, photon_energy, alpha, unit_wavevector, omega)

Create a [`LaserBeam`](@ref) object from a set of physically meaningful parameters.

This high-level constructor allows flexible initialization from either the laser’s
field strength (`a0`), its intensity, or any of the equivalent photon energy representations
(`wavelength`, `frequency`, or `photon_energy`).
If some parameters are missing, they are inferred where possible.

# Keyword Arguments

* `a0::Real = nothing`: Classical nonlinearity parameter (dimensionless).
  If not provided but `intensity` and `omega` (or its equivalents) are given, it is inferred.
* `intensity::Quantity = nothing`: Laser intensity (e.g. `1e20u"W/cm^2"`).
  Used to infer `a0` if not explicitly specified.
* `wavelength::Quantity = nothing`: Laser wavelength (e.g. `800u"nm"`).
  Used to infer photon energy if `frequency` or `photon_energy` are not given.
* `frequency::Quantity = nothing`: Laser frequency (alternative to wavelength).
  Used to infer photon energy if `wavelength` or `photon_energy` are not given.
* `photon_energy::Quantity = nothing`: Photon energy in eV (preferred).
  Used directly if provided; otherwise inferred from `frequency` or `wavelength`.
* `omega::Quantity = photon_energy`: Alias for `photon_energy`, for consistency with natural-unit notation.
* `alpha::Real = nothing`: Polar angle of the laser propagation direction (in radians).
  If not provided, it is derived from `unit_wavevector` if available.
* `unit_wavevector::AbstractVector = nothing`: 3D propagation direction of the laser.
  Used to infer the polar angle `alpha`.

# Behavior

At least one of `photon_energy`, `frequency`, or `wavelength` must be specified.
If `a0` is omitted but `intensity` is provided, it is inferred.

# Returns

A fully constructed [`LaserBeam`](@ref) object with:

* and a valid ``a_0`` parameter.
* propagation polar angle (``\\alpha``),
* photon energy ``\\omega`` (converted to eV),

# Throws

* `ArgumentError` if insufficient information is provided to determine photon energy.

# Example

```julia
# Construct from wavelength and intensity
L1 = laserbeam(wavelength = 800u"nm", intensity = 1e20u"W/cm^2", alpha = π/3)

# Construct from photon energy and direction vector
L2 = laserbeam(photon_energy = 1.55u"eV", unit_wavevector = [0, 0, 1], a0 = 1.0)

# Construct from frequency and intensity
L3 = laserbeam(frequency = 3.75e14u"Hz", intensity = 5e19u"W/cm^2")
```

"""
function laserbeam(;
        a0 = nothing,
        intensity = nothing,
        wavelength = nothing,
        frequency = nothing, #
        photon_energy = nothing,
        alpha = nothing,
        unit_wavevector = nothing,
        omega = photon_energy,
    )

    if alpha == nothing && unit_wavevector != nothing
        k = normalize(SVector{3}(unit_wavevector))
        alpha = _safe_polar_angle(k)
    end

    # infer photon energy
    if omega == nothing
        if wavelength == nothing && frequency != nothing
            omega = mynatural(Unitful.h * frequency)
        elseif frequency == nothing && wavelength != nothing
            omega = mynatural(Unitful.h * Unitful.c0 / wavelength)
        else
            throw(
                ArgumentError(
                    "laser beams must be initialized with one of either photon_energy, frequency, or wavelength"
                )
            )
        end
    else
        omega = omega |> u"eV"
    end


    if a0 == nothing && intensity != nothing
        a0 = _compute_a0_from_omega(intensity, omega)
    end

    return LaserBeam(a0, omega, alpha)
end

### accessors
"""
    photon_energy(L::LaserBeam)

Return the photon energy in eV associated with the laser frequency.
"""
photon_energy(L::LaserBeam) = mynatural(frequency(L))

"""
    classical_nonlinearity_parameter(L::LaserBeam)

Return the dimensionless laser strength parameter ``a_0``.
"""
classical_nonlinearity_parameter(L::LaserBeam) = L.a0

"""
    polar_angle(L::LaserBeam)

Return polar angle w.r.t. the z-axis of the laser propagation direction.
"""
polar_angle(L::LaserBeam) = L.alpha

### Derived quantities
#
"""
    unit_wavevector(L::LaserBeam)

Return the normalized propagation direction of the laser beam as a 3D vector.
"""
unit_wavevector(L::LaserBeam) = _unit_vec(L.alpha)

"""
    wavelength(L::LaserBeam)

Return the wavelength associated with the laser beam (in ``\\mu\\mathrm{m}``}).
"""
function wavelength(L::LaserBeam)
    photon_energy = uconvert(u"J", L.omega)
    lambda = Unitful.c0 * Unitful.h / photon_energy |> u"µm"
    return lambda
end

"""
    intensity(L::LaserBeam)

Compute the laser intensity from `a0` and `wavelength`.
Returned as a physical quantity in ``\\mathrm{W}/\\mathrm{cm}^2``.
"""
intensity(L::LaserBeam) = _compute_intensity_from_omega(L.a0, L.omega)

"""
    frequency(L::LaserBeam)

Return the frequency of the laser in Hz, defined as ``\\nu = c/\\lambda``.

Note: This uses the linear frequency (Hz), not the angular frequency (rad/s).
"""
frequency(L::LaserBeam, unit::Unitful.Units = u"Hz") = uconvert(u"J", L.omega) / Unitful.h |> unit

"""
    four_momentum(L::LaserBeam)

Return the photon four-momentum associated with the laser beam:

``k^\\mu = (\\omega, \\vec{k})``

where ``\\omega`` is the photon energy in natural units and ``\\vec{k}`` points along the normalized propagation direction.
"""
function four_momentum(L::LaserBeam{T}) where {T}
    omega = photon_energy(L)
    unit_three_mom = unit_wavevector(L)
    return omega * SLorentzVector(one(T), unit_three_mom...)
end
