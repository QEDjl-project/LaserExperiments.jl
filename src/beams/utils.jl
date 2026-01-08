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
