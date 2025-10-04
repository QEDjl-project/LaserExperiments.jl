using Test
using Unitful
using Random
using QEDbase
using QEDcore
using StaticArrays
using LinearAlgebra
using LaserExperiments  # your package

RNG = Xoshiro(161)

A0_ARR = (0.0, 1.0e-4, 1.0, 10.0, 100.0)
ALPHA_ARR = (0.0, rand(RNG) * pi, Float64(pi))
OMEGA_ARR = (1.0e-2, 1.0, rand(RNG), 1.0e3, 1.0e6) .* 1u"eV"

_frequency_from_energy(omega) = omega / Unitful.h |> u"Hz"
_wavelength_from_energy(omega) = Unitful.h * Unitful.c0 / omega |> u"µm"
function _intensity_from_a0(a0, omega)
    omega_over_eV = ustrip(uconvert(u"eV", omega))
    return (a0 * omega_over_eV / 7.5)^2 * 1.0e20u"W/cm^2"
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

@inline function _check_laser_beam(laser, a0, omega, alpha)
    @test isapprox(classical_nonlinearity_parameter(laser), a0)
    @test unit(classical_nonlinearity_parameter(laser)) == NoUnits

    @test isapprox(photon_energy(laser), omega)
    @test unit(photon_energy(laser)) == u"eV"

    @test isapprox(polar_angle(laser), alpha)
    return @test unit(polar_angle(laser)) == NoUnits
end

@testset "a0 = $A0, omega=$OMEGA, alpha=$ALPHA" for (A0, OMEGA, ALPHA) in Iterators.product(
        A0_ARR, OMEGA_ARR, ALPHA_ARR
    )

    @testset "LaserBeam basic construction" begin
        # Minimal construction with explicit parameters
        L = LaserBeam(A0, OMEGA, ALPHA)
        _check_laser_beam(L, A0, OMEGA, ALPHA)

        # Check unit of omega
        L2 = LaserBeam(A0, uconvert(u"J", OMEGA), ALPHA)
        _check_laser_beam(L2, A0, OMEGA, ALPHA)
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
        @test isapprox(four_momentum(L), FOUR_MOMENTUM)
        #@test all(unit.(four_momentum(L)) .== u"eV")
    end
end

@testset "invalid input" begin

    # wrong a0
    @test_throws ArgumentError LaserBeam(-1.0, 1.0u"eV", 0.0)

    # wrong omega
    @test_throws ArgumentError LaserBeam(1.0, -1.0u"eV", 0.0)

    # wrong alpha
    @test_throws ArgumentError LaserBeam(1.0, 1.0u"eV", -1.0)
    @test_throws ArgumentError LaserBeam(1.0, 1.0u"eV", 1.1 * pi)
end
