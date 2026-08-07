function _check_isapprox_quantity(a, b, kwargs...)
    return if ismissing(b)
        @test ismissing(a)
    else
        @test isapprox(ustrip(a), ustrip(b), kwargs...)
        @test unit(a) == unit(b)
    end
end

function _equal_or_missing(a, b)
    if ismissing(b)
        return ismissing(a)
    else
        return ismissing(a) ? false : a == b
    end
end

@inline function _check_laser_beam(laser_params, a0, ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio)
    _check_isapprox_quantity(classical_nonlinearity_parameter(laser_params), a0)
    _check_isapprox_quantity(photon_energy(laser_params), ph_en)
    _check_isapprox_quantity(pulse_energy(laser_params), pu_en)
    _check_isapprox_quantity(pulse_duration(laser_params), pu_du)
    _check_isapprox_quantity(repetition_rate(laser_params), rep_rate)
    _check_isapprox_quantity(coherence_length(laser_params), coh_len)
    @test _equal_or_missing(polarization_ratio(laser_params), pol_ratio)

    return nothing
end


_frequency_from_energy(omega) = omega / Unitful.h |> u"Hz"
_wavelength_from_energy(omega) = Unitful.h * Unitful.c0 / omega |> u"µm"
function _intensity_from_a0(a0, omega)
    omega_over_eV = ustrip(uconvert(u"eV", omega))
    return (a0 * omega_over_eV / 7.5)^2 * 1.0e20u"W/cm^2"
end
_average_power_from_pulse_energy(pulse_energy, repetition_rate) = pulse_energy * repetition_rate
