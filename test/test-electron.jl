using Test
using Unitful
using Random
using QEDcore
using StaticArrays
using LinearAlgebra
using LaserExperiments

RNG = Xoshiro(161)

ALPHA_ARR = (0.0, rand(RNG) * pi, Float64(pi))
ENERGY_ARR = ME_eV .+ (1.0e-2, 1.0, rand(RNG), 1.0e3, 1.0e6) .* 1u"eV" .|> u"eV"

@inline function _check_electron_beam(elec_beam, energy, alpha)
    @test isapprox(total_energy(elec_beam), energy)
    @test unit(total_energy(elec_beam)) == u"eV"

    @test isapprox(polar_angle(elec_beam), alpha)
    @test unit(polar_angle(elec_beam)) == NoUnits

    return nothing
end

@testset "energy=$ENERGY, alpha=$ALPHA" for (ENERGY, ALPHA) in Iterators.product(
        ENERGY_ARR, ALPHA_ARR
    )

    @testset "ElectronBeam basic construction" begin
        # Minimal construction with explicit parameters
        ELEC_BEAM = ElectronBeam(ENERGY, ALPHA)
        _check_electron_beam(ELEC_BEAM, ENERGY, ALPHA)

        # Check unit of omega
        ELEC_BEAM2 = ElectronBeam(uconvert(u"J", ENERGY), ALPHA)
        _check_electron_beam(ELEC_BEAM2, ENERGY, ALPHA)
    end

    @testset "ElectronBeam keyword constructor inference" begin

        @testset "energy" begin
            UNIT_DIR = _unit_vec(ALPHA)

            b1 = electronbeam(energy = ENERGY, alpha = ALPHA)
            _check_electron_beam(b1, ENERGY, ALPHA)

            b2 = electronbeam(energy = ENERGY, unit_direction = UNIT_DIR)
            _check_electron_beam(b2, ENERGY, ALPHA)
        end

        @testset "Lorentz factor" begin
            GAMMA = ustrip(ENERGY / ME_eV)
            UNIT_DIR = _unit_vec(ALPHA)

            b1 = electronbeam(gamma = GAMMA, alpha = ALPHA)
            _check_electron_beam(b1, _energy_from_gamma(GAMMA), ALPHA)

            b2 = electronbeam(gamma = GAMMA, unit_direction = UNIT_DIR)
            _check_electron_beam(b2, _energy_from_gamma(GAMMA), ALPHA)
        end

        @testset "beta" begin
            BETA = ustrip(sqrt(1 - (ME_eV / ENERGY)^2))
            UNIT_DIR = _unit_vec(ALPHA)

            b1 = electronbeam(beta = BETA, alpha = ALPHA)
            _check_electron_beam(b1, _energy_from_beta(BETA), ALPHA)

            b2 = electronbeam(beta = BETA, unit_direction = UNIT_DIR)
            _check_electron_beam(b2, _energy_from_beta(BETA), ALPHA)
        end

        @testset "momentum" begin
            P_MAG = sqrt(ENERGY^2 - ME_eV^2)
            UNIT_DIR = _unit_vec(ALPHA)

            b1 = electronbeam(p_mag = P_MAG, alpha = ALPHA)
            _check_electron_beam(b1, _energy_from_pmag(P_MAG), ALPHA)

            b2 = electronbeam(p_mag = P_MAG, unit_direction = UNIT_DIR)
            _check_electron_beam(b2, _energy_from_pmag(P_MAG), ALPHA)
        end
    end

    @testset "Derived quantities" begin
        B = ElectronBeam(ENERGY, ALPHA)

        # Groundtruths
        GAMMA = ustrip(ENERGY / ME_eV)
        BETA = ustrip(sqrt((1 - ME_eV / ENERGY) * (1 + ME_eV / ENERGY)))
        P_MAG = sqrt((ENERGY - ME_eV) * (ENERGY + ME_eV)) |> u"eV"
        UNIT_DIR = _unit_vec(ALPHA)
        FOUR_MOME_eVNTUM = SLorentzVector(ENERGY, P_MAG * UNIT_DIR...)

        @test isapprox(lorentz_factor(B), GAMMA)
        @test isapprox(beta(B), BETA)
        @test isapprox(momentum_magnitude(B), P_MAG)
        @test isapprox(unit_direction(B), UNIT_DIR)
        @test isapprox(four_momentum(B), FOUR_MOME_eVNTUM)
    end
end

@testset "invalid input" begin
    # wrong energy
    @test_throws ArgumentError ElectronBeam(-1.0u"eV", 0.0)

    # wrong alpha
    @test_throws ArgumentError ElectronBeam(1.0u"MeV", -0.1)
    @test_throws ArgumentError ElectronBeam(1.0u"MeV", 3.5)
end
