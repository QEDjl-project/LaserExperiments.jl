abstract type AbstractElectronBeam <: AbstractBeam end

"""
    struct ElectronBeam{T} <: AbstractElectronBeam
    ElectronBeam(energy, alpha)

Representation of a relativistic electron beam characterized by its total energy
and polar propagation angle.

## Fields
- `energy::Quantity{T}`: Total electron energy (in eV).
- `alpha::T`: Polar angle (in radians) defining the propagation direction.

## Construction

A new `ElectronBeam` can be created via the outer constructor [`electronbeam`](@ref).
The inner constructor performs sanity checks and unit normalization.

## Example

```

E = ElectronBeam(10.0u"GeV", π/6)

```
"""
struct ElectronBeam{T} <: AbstractElectronBeam
    energy::Quantity{T}
    alpha::T

    function ElectronBeam(energy::Quantity{T}, alpha::T) where {T <: Real}
        natural_energy = mynatural(energy)
        ustrip(natural_energy) >= T(ustrip(ME_eV))|| throw(ArgumentError("electron energy must be larger than the electron mass"))
        zero(T) <= alpha <= T(pi) || throw(ArgumentError("polar angle `alpha` must be between 0 and pi"))
        return new{T}(natural_energy, alpha)
    end
end


### Extended Constructor

"""
    electronbeam(; energy, gamma, beta, p_mag, alpha, unit_direction)

Create an [`ElectronBeam`](@ref) object from physically meaningful parameters.

This high-level constructor allows flexible initialization of an electron beam
from its total energy, Lorentz factor (`gamma`), velocity (`beta`), or momentum magnitude (`p_mag`).
If some parameters are missing, they are inferred self-consistently from the relativistic on-shell relations.

# Keyword Arguments

* `energy::Quantity = nothing`: Total energy (including rest mass), e.g. `10u"GeV"`.
* `gamma::Real = nothing`: Lorentz factor, ``\\gamma = E / (ME_eV c^2)``.
* `beta::Real = nothing`: Dimensionless velocity, ``\\beta = \\sqrt{1 - 1/\\gamma^2}``.
* `p_mag::Quantity = nothing`: Magnitude of three-momentum (in eV, natural units).
* `alpha::Real = nothing`: Polar angle (radians) defining propagation direction.
  If not given, it is inferred from `unit_direction`.
* `unit_direction::AbstractVector = nothing`: 3D propagation direction vector.
  Used to infer `alpha`.

# Behavior

At least one of `energy`, `gamma`, or `p_mag` must be specified.
If several are given, they are cross-checked for physical consistency.

# Returns

A fully constructed [`ElectronBeam`](@ref) with:
- total energy ``E``,
- polar angle ``\\alpha``.

# Throws

* `ArgumentError` if insufficient information is provided to determine energy.
* `ArgumentError` if provided parameters are inconsistent.

# Example

```julia
# Construct from total energy and polar angle
E1 = electronbeam(energy = 10u"GeV", alpha = π/3)

# Construct from gamma factor
E2 = electronbeam(gamma = 2e4, unit_direction = [0, 0, 1])

# Construct from momentum magnitude
E3 = electronbeam(p_mag = 10u"GeV", alpha = π/6)
```

"""
function electronbeam(;
        energy = nothing,
        gamma = nothing,
        beta = nothing,
        p_mag = nothing,
        alpha = nothing,
        unit_direction = nothing
    )
    # Infer propagation direction
    if alpha == nothing
        if unit_direction == nothing
            throw(ArgumentError("either `alpha` or `unit_direction` must be provided"))
        end
        k = normalize(SVector{3}(unit_direction))
        alpha = _safe_polar_angle(k)
    end

    # Ensure at least one kinematic quantity is provided
    if all(x -> x === nothing, (energy, gamma, beta, p_mag))
        throw(ArgumentError("must provide one of `energy`, `gamma`, `beta`, or `p_mag`"))
    end


    if energy == nothing
        if gamma !== nothing && beta == nothing && p_mag == nothing
            # --- Case 1: start from gamma ---
            energy = mynatural(gamma * ME_eV)

        elseif beta !== nothing && gamma == nothing && p_mag == nothing
            # --- Case 2: start from beta ---
            gamma = 1 / sqrt((1 - beta) * (1 + beta))
            energy = mynatural(gamma * ME_eV)

        elseif p_mag !== nothing && beta == nothing && gamma == nothing
            # --- Case 3: start from p_mag ---
            p_mag = ustrip(mynatural(p_mag))
            me = ustrip(ME_eV)
            #E_val = sqrt((ustrip(p_mag))^2 + (ustrip(ME_eV))^2) * u"eV"
            energy = sqrt(p_mag^2 + me^2) * u"eV"
        else
            throw(
                ArgumentError(
                    "provide exactly one of `energy`, `gamma`, `beta`, or `p_mag`"
                )
            )
        end
    end

    return ElectronBeam(energy, alpha)

end

### Accessors

"""
    total_energy(E::ElectronBeam)

Return the total energy (including rest mass) of the electron beam in eV.
"""
total_energy(E::ElectronBeam) = E.energy

"""
    polar_angle(E::ElectronBeam)

Return the polar angle (radians) of the electron propagation direction.
"""
polar_angle(E::ElectronBeam) = E.alpha

### Derived quantities

"""
    lorentz_factor(E::ElectronBeam)

Return the Lorentz factor `\\gamma = E / (ME_eV c^2)`.
"""
lorentz_factor(E::ElectronBeam) = ustrip(total_energy(E) / ME_eV)

"""
    beta(E::ElectronBeam)

Return the dimensionless velocity `\\beta = \\sqrt{1 - 1/\\gamma^2}`.
"""
function beta(E::ElectronBeam)
    inv_gamma = inv(lorentz_factor(E))
    return sqrt((1 - inv_gamma) * (1 + inv_gamma))
end

"""
    momentum_magnitude(E::ElectronBeam)

Return the absolute value of the three-momentum |p| in eV (natural units),
computed from the on-shell relation `p = \\sqrt{E^2 - ME_eV^2}`.
"""
function momentum_magnitude(E::ElectronBeam)

    #E_eV = ustrip(total_energy(E))
    En = ustrip(total_energy(E))
    me = ustrip(ME_eV)
    #ME_eV = ustrip(ME)
    return sqrt((En - me) * (En + me)) * u"eV"
end

"""
    unit_direction(E::ElectronBeam)

Return the normalized propagation direction of the electron beam as a 3D vector.
"""
unit_direction(E::ElectronBeam) = _unit_vec(E.alpha)

"""
    four_momentum(E::ElectronBeam)

Return the four-momentum of the electron beam in natural units:

`p^\\mu = (E, \\vec{p})`

where `E` is the total energy and `\\vec{p}` points along the propagation direction.
"""
function four_momentum(E::ElectronBeam{T}) where {T}
    E_tot = total_energy(E)
    p_mag = momentum_magnitude(E)
    unit_vec = unit_direction(E)
    return SLorentzVector(E_tot, p_mag * unit_vec...)
end
