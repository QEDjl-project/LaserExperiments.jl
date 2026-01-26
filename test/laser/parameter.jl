using Test
using Unitful
using Random
using QEDcore
using StaticArrays
using LinearAlgebra
using LaserExperiments

include("utils.jl")

RNG = Xoshiro(161)

A0S = (rand(RNG),)
PHOTON_ENERGIES = (rand(RNG) * u"eV",)
PULSE_ENERGIES = (rand(RNG) * u"J", missing)
PULSE_DURATIONS = (rand(RNG) * u"fs", missing)
REP_RATES = (rand(RNG) * u"Hz", missing)
COHERENCE_LENGTHS = (rand(RNG) * u"cm", missing)
POL_RATIOS = [
    (1.0e3 * rand(RNG), 1.0),
    missing,
]


@testset "a0: $a0, ph_en: $ph_en, pu_en: $pu_en, pu_du: $pu_du, rep_rate: $rep_rate, coh_len: $coh_len, pol_ratio: $pol_ratio" for (a0, ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio) in Iterators.product(
        A0S, PHOTON_ENERGIES, PULSE_ENERGIES, PULSE_DURATIONS, REP_RATES, COHERENCE_LENGTHS, POL_RATIOS
    )

    @testset "basic construction" begin
        # Minimal construction with explicit parameters
        L = LaserParameters(a0, ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio)
        _check_laser_beam(L, a0, ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio)
    end

    @testset "invalid input" begin

        # wrong a0
        if !ismissing(a0)
            @test_throws ArgumentError LaserParameters(-a0, ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio)
            @test_throws ArgumentError LaserParameters(zero(a0), ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio)
        end

        # wrong photon_energy
        if !ismissing(ph_en)
            @test_throws ArgumentError LaserParameters(a0, -ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio)
            @test_throws ArgumentError LaserParameters(a0, zero(ph_en), pu_en, pu_du, rep_rate, coh_len, pol_ratio)
        end

        # wrong pulse_energy
        if !ismissing(pu_en)
            @test_throws ArgumentError LaserParameters(a0, ph_en, -pu_en, pu_du, rep_rate, coh_len, pol_ratio)
            @test_throws ArgumentError LaserParameters(a0, ph_en, zero(pu_en), pu_du, rep_rate, coh_len, pol_ratio)
        end

        # wrong pulse_duration
        if !ismissing(pu_du)
            @test_throws ArgumentError LaserParameters(a0, ph_en, pu_en, -pu_du, rep_rate, coh_len, pol_ratio)
            @test_throws ArgumentError LaserParameters(a0, ph_en, pu_en, zero(pu_du), rep_rate, coh_len, pol_ratio)
        end

        # wrong rep_rate
        if !ismissing(rep_rate)
            @test_throws ArgumentError LaserParameters(a0, ph_en, pu_en, pu_du, -rep_rate, coh_len, pol_ratio)
            @test_throws ArgumentError LaserParameters(a0, ph_en, pu_en, pu_du, zero(rep_rate), coh_len, pol_ratio)
        end

        # wrong coh_len
        if !ismissing(coh_len)
            @test_throws ArgumentError LaserParameters(a0, ph_en, pu_en, pu_du, rep_rate, -coh_len, pol_ratio)
            @test_throws ArgumentError LaserParameters(a0, ph_en, pu_en, pu_du, rep_rate, zero(coh_len), pol_ratio)
        end

        # wrong pol_ratio
        if !ismissing(pol_ratio)
            @test_throws ArgumentError LaserParameters(a0, ph_en, pu_en, pu_du, rep_rate, coh_len, (zero(pol_ratio[1]), pol_ratio[2]))
            @test_throws ArgumentError LaserParameters(a0, ph_en, pu_en, pu_du, rep_rate, coh_len, (pol_ratio[1], zero(pol_ratio[2])))
        end
    end

end

@testset "LaserParameters keyword constructor" begin
    @testset "defaults" begin
        @testset "a0: $a0, ph_en: ph_en" for (a0, ph_en) in Iterators.product(A0S, PHOTON_ENERGIES)
            L = laserparameter(a0 = a0, photon_energy = ph_en)
            _check_laser_beam(L, a0, ph_en, missing, missing, missing, missing, missing)
        end
    end

    @testset "a0: $a0, ph_en: $ph_en, pu_en: $pu_en, pu_du: $pu_du, rep_rate: $rep_rate, coh_len: $coh_len, pol_ratio: $pol_ratio" for (a0, ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio) in Iterators.product(
            A0S, PHOTON_ENERGIES, PULSE_ENERGIES, PULSE_DURATIONS, REP_RATES, COHERENCE_LENGTHS, POL_RATIOS
        )
        @testset "intensity" begin
            L = laserparameter(;
                a0 = a0,
                photon_energy = ph_en,
                pulse_energy = pu_en,
                pulse_duration = pu_du,
                repetition_rate = rep_rate,
                coherence_length = coh_len,
                polarization_ratio = pol_ratio,
            )
            _check_laser_beam(L, a0, ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio)

            INTENSITY = _intensity_from_a0(a0, ph_en)
            L = laserparameter(;
                intensity = INTENSITY,
                photon_energy = ph_en,
                pulse_energy = pu_en,
                pulse_duration = pu_du,
                repetition_rate = rep_rate,
                coherence_length = coh_len,
                polarization_ratio = pol_ratio,
            )
            _check_laser_beam(L, a0, ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio)
        end

        @testset "photon energy" begin
            FREQUENCY = _frequency_from_energy(ph_en)
            L = laserparameter(;
                a0 = a0,
                frequency = FREQUENCY,
                pulse_energy = pu_en,
                pulse_duration = pu_du,
                repetition_rate = rep_rate,
                coherence_length = coh_len,
                polarization_ratio = pol_ratio,
            )
            _check_laser_beam(L, a0, ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio)

            WAVELENGTH = _wavelength_from_energy(ph_en)
            L = laserparameter(;
                a0 = a0,
                wavelength = WAVELENGTH,
                pulse_energy = pu_en,
                pulse_duration = pu_du,
                repetition_rate = rep_rate,
                coherence_length = coh_len,
                polarization_ratio = pol_ratio,
            )
            _check_laser_beam(L, a0, ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio)
        end
        @testset "pulse energy" begin
            if (ismissing(pu_en) && ismissing(rep_rate)) || (!ismissing(pu_en) && !ismissing(rep_rate))
                AVERAGE_POWER = _average_power_from_pulse_energy(pu_en, rep_rate)
                L = laserparameter(;
                    a0 = a0,
                    omega = ph_en,
                    average_power = AVERAGE_POWER,
                    pulse_duration = pu_du,
                    repetition_rate = rep_rate,
                    coherence_length = coh_len,
                    polarization_ratio = pol_ratio,
                )
                _check_laser_beam(L, a0, ph_en, pu_en, pu_du, rep_rate, coh_len, pol_ratio)
            end
        end

        @testset "derived quantities" begin
            L = laserparameter(;
                a0 = a0,
                photon_energy = ph_en,
                pulse_energy = pu_en,
                pulse_duration = pu_du,
                repetition_rate = rep_rate,
                coherence_length = coh_len,
                polarization_ratio = pol_ratio,
            )

            INTENSITY = _intensity_from_a0(a0, ph_en)
            @test isapprox(intensity(L), INTENSITY)
            FREQUENCY = _frequency_from_energy(ph_en)
            @test isapprox(frequency(L), FREQUENCY)
            WAVELENGTH = _wavelength_from_energy(ph_en)
            @test isapprox(wavelength(L), WAVELENGTH)
            if (ismissing(pu_en) && ismissing(rep_rate)) || (!ismissing(pu_en) && !ismissing(rep_rate))
                AVERAGE_POWER = _average_power_from_pulse_energy(pu_en, rep_rate)
                _check_isapprox_quantity(average_power(L), AVERAGE_POWER)
            end
        end
    end
end
