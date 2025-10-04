function _unit_vec(polar, azimuth = zero(polar))
    sth, cth = sincos(polar)
    sphi, cphi = sincos(azimuth)
    return SVector(
        cphi * sth,
        sphi * sth,
        cth
    )
end

_frequency_from_energy(omega) = omega / Unitful.h |> u"Hz"
_wavelength_from_energy(omega) = Unitful.h * Unitful.c0 / omega |> u"µm"
function _intensity_from_a0(a0, omega)
    omega_over_eV = ustrip(uconvert(u"eV", omega))
    return (a0 * omega_over_eV / 7.5)^2 * 1.0e20u"W/cm^2"
end

_energy_from_pmag(p_mag) = sqrt(p_mag^2 + ME_eV^2) |> u"eV"
_energy_from_gamma(gamma) = gamma * ME_eV |> u"eV"
_energy_from_beta(beta) = ME_eV / sqrt((1 - beta) * (1 + beta)) |> u"eV"
