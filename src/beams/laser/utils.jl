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
