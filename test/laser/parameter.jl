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

    #=
        laserparameter(;
                          a0 = nothing,
                          intensity = nothing,

                          wavelength = nothing,
                          frequency = nothing,
                          photon_energy = nothing,
                          omega = photon_energy,

                          pulse_duration = nothing,
                          coherence_length = nothing,
                          polarization_ratio = (one(Float64), one(Float64)),

                          pulse_energy = nothing,
                          average_power = nothing,
                          repetition_rate = nothing,
                          )

    end
    @testset "LaserBeam keyword constructor inference" begin

        @testset "photon_energy" begin
            INTENSITY = _intensity_from_a0(A0, OMEGA)
            UNIT_WAVEVECTOR = _unit_vec(ALPHA)

            L = laserbeam(a0 = A0, photon_energy = OMEGA, alpha = ALPHA)
            _check_laser_beam(L, A0, OMEGA, ALPHA)

            L2 = laserbeam(intensity = INTENSITY, photon_energy = OMEGA, alpha = ALPHA)
            _check_laser_beam(L2, A0, OMEGA, ALPHA)

            L3 = laserbeam(a0 = A0, photon_energy = OMEGA, unit_wavevector = UNIT_WAVEVECTOR)
            _check_laser_beam(L3, A0, OMEGA, ALPHA)

            L4 = laserbeam(intensity = INTENSITY, photon_energy = OMEGA, unit_wavevector = UNIT_WAVEVECTOR)
            _check_laser_beam(L4, A0, OMEGA, ALPHA)

            # alias: photon_energy==omega
            L5 = laserbeam(a0 = A0, omega = OMEGA, alpha = ALPHA)
            _check_laser_beam(L5, A0, OMEGA, ALPHA)

            L6 = laserbeam(intensity = INTENSITY, omega = OMEGA, alpha = ALPHA)
            _check_laser_beam(L6, A0, OMEGA, ALPHA)

            L7 = laserbeam(a0 = A0, omega = OMEGA, unit_wavevector = UNIT_WAVEVECTOR)
            _check_laser_beam(L7, A0, OMEGA, ALPHA)

            L8 = laserbeam(intensity = INTENSITY, omega = OMEGA, unit_wavevector = UNIT_WAVEVECTOR)
            _check_laser_beam(L8, A0, OMEGA, ALPHA)
        end

        @testset "frequency" begin
            FREQUENCY = _frequency_from_energy(OMEGA)
            INTENSITY = _intensity_from_a0(A0, OMEGA)
            UNIT_WAVEVECTOR = _unit_vec(ALPHA)

            L = laserbeam(a0 = A0, frequency = FREQUENCY, alpha = ALPHA)
            _check_laser_beam(L, A0, OMEGA, ALPHA)

            L2 = laserbeam(intensity = INTENSITY, frequency = FREQUENCY, alpha = ALPHA)
            _check_laser_beam(L2, A0, OMEGA, ALPHA)

            L3 = laserbeam(a0 = A0, frequency = FREQUENCY, unit_wavevector = UNIT_WAVEVECTOR)
            _check_laser_beam(L3, A0, OMEGA, ALPHA)

            L4 = laserbeam(intensity = INTENSITY, frequency = FREQUENCY, unit_wavevector = UNIT_WAVEVECTOR)
            _check_laser_beam(L4, A0, OMEGA, ALPHA)
        end

        @testset "wavelength" begin
            WAVELENGTH = _wavelength_from_energy(OMEGA)
            INTENSITY = _intensity_from_a0(A0, OMEGA)
            UNIT_WAVEVECTOR = _unit_vec(ALPHA)

            L = laserbeam(a0 = A0, wavelength = WAVELENGTH, alpha = ALPHA)
            _check_laser_beam(L, A0, OMEGA, ALPHA)

            L2 = laserbeam(intensity = INTENSITY, wavelength = WAVELENGTH, alpha = ALPHA)
            _check_laser_beam(L2, A0, OMEGA, ALPHA)

            L3 = laserbeam(a0 = A0, wavelength = WAVELENGTH, unit_wavevector = UNIT_WAVEVECTOR)
            _check_laser_beam(L3, A0, OMEGA, ALPHA)

            L4 = laserbeam(intensity = INTENSITY, wavelength = WAVELENGTH, unit_wavevector = UNIT_WAVEVECTOR)
            _check_laser_beam(L4, A0, OMEGA, ALPHA)
        end
    end

    @testset "Derived quantities" begin
        L = LaserBeam(A0, OMEGA, ALPHA)

        # groundtruths
        FREQUENCY = _frequency_from_energy(OMEGA)
        WAVELENGTH = _wavelength_from_energy(OMEGA)
        INTENSITY = _intensity_from_a0(A0, OMEGA)
        UNIT_WAVEVECTOR = _unit_vec(ALPHA)
        FOUR_MOMENTUM = OMEGA * SLorentzVector(1.0, _unit_vec(ALPHA)...)

        @test isapprox(frequency(L), FREQUENCY)
        @test unit(frequency(L)) == u"Hz"

        @test isapprox(wavelength(L), WAVELENGTH)
        @test unit(wavelength(L)) == u"µm"

        @test isapprox(intensity(L), INTENSITY)
        @test unit(intensity(L)) == u"W/cm^2"

        @test isapprox(unit_wavevector(L), UNIT_WAVEVECTOR)

        # FIXME: implement extension for QEDcore to use Unitful.Quantity properly
        #@test isapprox(four_momentum(L), FOUR_MOMENTUM)
        #@test all(unit.(four_momentum(L)) .== u"eV")
    end
    =#
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
    end
end
#=
@testset "invalid input" begin

    # wrong a0
    @test_throws ArgumentError LaserBeam(-1.0, 1.0u"eV", 0.0)

    # wrong omega
    @test_throws ArgumentError LaserBeam(1.0, -1.0u"eV", 0.0)

    # wrong alpha
    @test_throws ArgumentError LaserBeam(1.0, 1.0u"eV", -1.0)
    @test_throws ArgumentError LaserBeam(1.0, 1.0u"eV", 1.1 * pi)
end
=#
