# TODO: use internal unit system to avoid these constants

abstract type AbstractLaserBeam <: AbstractBeam end

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

function _unit_vec(polar, azimuth = zero(polar))
    sth, cth = sincos(polar)
    sphi, cphi = sincos(azimuth)
    return SVector(
        cphi * sth,
        sphi * sth,
        cth
    )
end

function _safe_polar_angle(unit_vector::SVector{3, T}) where {T}

    x, y, z = unit_vector

    if iszero(z)
        return T(pi / 2)
    end

    rad = atan(sqrt(x^2 + y^2) / z)

    return z > zero(z) ? T(rad) : T(pi + rad)
end
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

# DiPiazza 2012
function _compute_a0_from_omega(intensity::Quantity{T}, omega::Quantity{T}) where {T}
    I = ustrip(uconvert(u"W/cm^2", intensity)) / 1.0e20
    omega_over_eV = ustrip(uconvert(u"eV", omega))
    return T(7.5) * sqrt(I) / omega_over_eV
end

function _compute_intensity_from_omega(a0::T, omega::Quantity{T}) where {T}
    omega_over_eV = ustrip(uconvert(u"eV", omega))
    return T((a0 * omega_over_eV / 7.5)^2 * 1.0e20) * 1u"W/cm^2"
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

polar_angle(L::LaserBeam) = L.alpha

### Derived quantities
#
"""
    unit_wavevector(L::LaserBeam)

Return the normalized propagation direction of the laser beam as a 3D vector.
"""
unit_wavevector(L::LaserBeam) = _unit_vec(L.alpha)

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
