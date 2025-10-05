function _unit_vec(polar, azimuth = zero(polar))
    sth, cth = sincos(polar)
    sphi, cphi = sincos(azimuth)
    return SVector(
        cphi * sth,
        sphi * sth,
        cth
    )
end

# see https://en.wikipedia.org/wiki/Spherical_coordinates
function _safe_polar_angle(unit_vector::SVector{3, T}) where {T}

    x, y, z = unit_vector

    if iszero(z)
        return T(pi / 2)
    end

    rad = atan(sqrt(x^2 + y^2) / z)

    return z > zero(z) ? T(rad) : T(pi + rad)
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
