using Test
using Unitful
using Random
using QEDbase
using QEDcore
using StaticArrays
using LinearAlgebra
using LaserExperiments  # your package

RNG = Xoshiro(161)

@testset "LaserBeam basic construction" begin
    # Minimal construction with explicit parameters
    A0 = rand(RNG) * rand(RNG, 1:100)
    WAVELENGTH = uconvert(u"m", rand(RNG) * rand(RNG, 1:10) * 1u"µm")
    UNIT_WAVEVECTOR = normalize(@SVector rand(RNG, 3))
    L = LaserBeam(A0, WAVELENGTH, UNIT_WAVEVECTOR)

    @test isapprox(classical_nonlinearity_parameter(L), A0)
    @test unit(classical_nonlinearity_parameter(L)) == NoUnits

    @test isapprox(wavelength(L), WAVELENGTH)
    @test unit(wavelength(L)) == u"µm"

    @test isapprox(unit_wavevector(L), UNIT_WAVEVECTOR)


    # Check normalization of direction
    L2 = laserbeam(a0 = A0, wavelength = WAVELENGTH, unit_wavevector = rand(RNG, 3))
    @test isapprox(norm(L2.unit_wavevector), 1.0)
end

@testset "LaserBeam keyword constructor inference" begin
    A0 = rand(RNG) * rand(RNG, 1:100)

    # wavelength from frequency
    FREQUENCY = rand(RNG) * 1.0e15u"Hz"
    L = laserbeam(a0 = A0, frequency = FREQUENCY)
    groundtruth = Unitful.c0 / FREQUENCY

    #@test abs(Unitful.c0 / FREQUENCY - L.wavelength) < 1e-10u"m"
    @test isapprox(groundtruth, L.wavelength)

    # wavelength from photon energy
    OMEGA = rand(RNG) * 10u"eV"
    L = laserbeam(a0 = A0, photon_energy = OMEGA)
    @test abs(Unitful.h * Unitful.c0 / OMEGA - L.wavelength) < 1.0e-10u"m"

    # a0 inferred from intensity and wavelength
    WAVELENGTH = rand(RNG) * rand(RNG, 1:10) * 1u"µm"
    INTENSITY = rand(RNG) * 1.0e24u"W/cm^2"
    L = laserbeam(intensity = INTENSITY, wavelength = WAVELENGTH)
    @test isapprox(L.a0, LaserExperiments.compute_a0(INTENSITY, WAVELENGTH))
end


@testset "Derived quantities" begin
    A0 = rand(RNG) * rand(RNG, 1:100)
    WAVELENGTH = rand(RNG) * rand(RNG, 1:10) * 1u"µm"
    UNIT_WAVEVECTOR = normalize(@SVector rand(RNG, 3))

    L = LaserBeam(A0, WAVELENGTH, UNIT_WAVEVECTOR)

    # Frequency = c / λ
    @test isapprox(frequency(L), Unitful.c0 / WAVELENGTH)

    # Photon energy = ħ * ω
    @test isapprox(photon_energy(L), LaserExperiments.mynatural(frequency(L)))

    # Intensity from a0 and λ
    @test isapprox(intensity(L), LaserExperiments.compute_intensity(A0, WAVELENGTH))

    # Four momentum checks
    p = four_momentum(L)

    # FIXME: implement extension for QEDcore to use Unitful.Quantity properly
    @test p isa QEDbase.AbstractLorentzVector
    @test isapprox(getE(p), photon_energy(L))
    @test isapprox(p[2:4], photon_energy(L) .* L.unit_wavevector)
end

@testset "compute_a0 and compute_intensity consistency" begin
    WAVELENGTH = rand(RNG) * rand(RNG, 1:10) * 1u"µm"
    INTENSITY = rand(RNG) * 1.0e24u"W/cm^2"
    a0 = LaserExperiments.compute_a0(INTENSITY, WAVELENGTH)
    I2 = LaserExperiments.compute_intensity(a0, WAVELENGTH)

    # Round-trip: intensity -> a0 -> intensity
    @test isapprox(I2, INTENSITY)
end
