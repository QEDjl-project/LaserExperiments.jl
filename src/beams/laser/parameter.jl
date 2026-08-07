abstract type AbstractLaserParameters end

struct LaserParameters{T, PhE, PuE, PuD, REP, COH, POL} <: AbstractLaserParameters
    a0::T
    photon_energy::PhE
    pulse_energy::PuE
    pulse_duration::PuD
    repetition_rate::REP
    coherence_length::COH
    polarization_ratio::POL  # e.g. (1, 100), maybe use irrationals

    function LaserParameters(
            a0::T,
            photon_energy::PhE,
            pulse_energy::PuE,
            pulse_duration::PuD,
            repetition_rate::REP,
            coherence_length::COH,
            polarization_ratio::POL,
        ) where {
            T <: Real,
            PhE <: Union{Unitful.Energy, Missing},
            PuE <: Union{Unitful.Energy, Missing},
            PuD <: Union{Unitful.Time, Missing},
            REP <: Union{Unitful.Frequency, Missing},
            COH <: Union{Unitful.Length, Missing},
            POL <: Union{NTuple{2, T}, Missing},
        }

        a0 > zero(T) || throw(ArgumentError("`a0` must be non-negative"))
        ismissing(photon_energy) || ustrip(photon_energy) > zero(T) || throw(ArgumentError("photon energy must be positive"))
        ismissing(pulse_energy) || ustrip(pulse_energy) > zero(T) || throw(ArgumentError("pulse energy must be positive"))
        ismissing(pulse_duration) || ustrip(pulse_duration) > zero(T) || throw(ArgumentError("pulse duration must be positive"))
        ismissing(repetition_rate) || ustrip(repetition_rate) > zero(T) || throw(ArgumentError("repetition rate must be positive"))
        ismissing(coherence_length) || ustrip(coherence_length) > zero(T) || throw(ArgumentError("coherence length must be positive"))
        ismissing(polarization_ratio) || all(p -> p > zero(T), polarization_ratio) ||
            throw(ArgumentError("polarization ratio entries must be positive"))

        return new{T, PhE, PuE, PuD, REP, COH, POL}(
            a0,
            photon_energy,
            pulse_energy,
            pulse_duration,
            repetition_rate,
            coherence_length,
            polarization_ratio
        )
    end
end

classical_nonlinearity_parameter(L::LaserParameters) = L.a0
photon_energy(L::LaserParameters) = L.photon_energy
pulse_energy(L::LaserParameters) = L.pulse_energy
pulse_duration(L::LaserParameters) = L.pulse_duration
repetition_rate(L::LaserParameters) = L.repetition_rate
coherence_length(L::LaserParameters) = L.coherence_length
polarization_ratio(L::LaserParameters) = L.polarization_ratio


function _infer_a0(a0, intensity, omega)
    if a0 === nothing
        intensity === nothing && throw(
            ArgumentError(
                "`a0` or `intensity` must be specified"
            )
        )
        a0 = _compute_a0_from_omega(intensity, omega)
    end
    return a0
end

function _infer_photon_energy(omega, wavelength, frequency)
    if omega === nothing
        if wavelength !== nothing
            omega = Unitful.h * Unitful.c0 / wavelength
        elseif frequency !== nothing
            omega = Unitful.h * frequency
        else
            throw(
                ArgumentError(
                    "`photon_energy`, `omega`, `wavelength` or `frequency` must be specified"
                )
            )

        end
    end
    return omega |> u"eV"
end

function _infer_pulse_energy(pulse_energy, average_power, repetition_rate)
    if pulse_energy === nothing
        average_power === nothing && (return missing)
        repetition_rate === nothing && throw(
            ArgumentError(
                "if `average_power` is given, `repetition_rate` must be specified as well"
            )
        )
        pulse_energy = average_power / repetition_rate
    end
    return ismissing(pulse_energy) ? pulse_energy : pulse_energy |> u"J"
end

function laserparameter(;
        a0 = nothing,
        intensity = nothing,

        wavelength = nothing,
        frequency = nothing,
        photon_energy = nothing,
        omega = photon_energy,

        pulse_duration = nothing,
        coherence_length = nothing,
        polarization_ratio = nothing,

        pulse_energy = nothing,
        average_power = nothing,
        repetition_rate = nothing,
    )

    # --- infer photon energy ---
    omega = _infer_photon_energy(omega, wavelength, frequency)

    # --- infer a0 ---
    a0 = _infer_a0(a0, intensity, omega)


    #pulse_energy = pulse_energy === nothing || ismissing(pulse_energy) ? missing : pulse_energy|> u"J"
    pulse_duration = pulse_duration === nothing || ismissing(pulse_duration) ? missing : pulse_duration |> u"fs"
    repetition_rate = repetition_rate === nothing || ismissing(repetition_rate) ? missing : repetition_rate |> u"Hz"
    coherence_length = coherence_length === nothing || ismissing(coherence_length) ? missing : coherence_length |> u"cm"

    # -- infer pulse_energy --
    pulse_energy = _infer_pulse_energy(pulse_energy, average_power, repetition_rate)

    polarization_ratio = isnothing(polarization_ratio) ? missing : polarization_ratio
    return LaserParameters(
        a0,
        omega,
        pulse_energy,
        pulse_duration,
        repetition_rate,
        coherence_length,
        polarization_ratio
    )
end

intensity(l::LaserParameters) = _compute_intensity_from_omega(l.a0, l.photon_energy)
frequency(l::LaserParameters) = l.photon_energy / Unitful.h |> u"Hz"
wavelength(l::LaserParameters) = Unitful.h * Unitful.c0 / l.photon_energy |> u"µm"
average_power(l::LaserParameters) = l.pulse_energy * l.repetition_rate
