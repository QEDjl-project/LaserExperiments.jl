# TODO:
# - refac this (make it more clear)
# - add some unit tests
#
const speed_of_light = 1u"c"
const hbar = Unitful.ħ
const hbarc = hbar * speed_of_light

function _get_dim_power(x::Unitful.Dimensions{T}) where {T}
    return prod(_get_dim_power.(T))
end

_get_dim_power(x::Unitful.Dimension{:Length}) = -x.power
_get_dim_power(x::Unitful.Dimension{:Time}) = -x.power
_get_dim_power(x::Unitful.Dimension{:Mass}) = x.power

function _get_dim_factor(x::Unitful.Dimensions{T}) where {T}
    return prod(@. _get_dim_factor(T, _get_dim_power(T)))
end

function _get_dim_factor(x::Unitful.Dimension{:Length}, p::Rational)
    return hbarc^p
end

function _get_dim_factor(x::Unitful.Dimension{:Mass}, p::Rational)
    return speed_of_light^(2 * p)
end

function _get_dim_factor(x::Unitful.Dimension{:Time}, p::Rational)
    return Unitful.h^(p)
end

function mynatural(q; base = u"eV")
    d = dimension(q)
    p = _get_dim_power(d)
    nunit = base^(p)
    nfac = _get_dim_factor(d)
    return uconvert(nunit, q * nfac)
end

# opt out, if the quanity is already in eV
gettypeparams(::Unitful.FreeUnits{T, U, V}) where {T, U, V} = T, U, V
const eV_DIM = gettypeparams(u"GeV")[2]

function mynatural(q::Quantity{T, eV_DIM}; base = u"eV") where {T <: Number}
    return uconvert(base, q)
end

# FIXME: generalize to powers `eV^p`, currently it errors for `p!=1

## conversion for structure factors

const me = ustrip(mynatural(Unitful.me))

function _eV2me_dimless(t_eV, p::Int = 1)
    return ustrip(t_eV) * me^p
end

function _energy2me_dimless(energy)
    energy_eV = mynatural(energy)
    return _eV2me_dimless(energy_eV, -1)
end

# simple unit transforms
# TODO:
# - write unit tests for this
# - if other unit transforms pop up, add them here
#
_me2eV(x_me, p = 1) = x_me * (ELECTRONMASS)^p
_me2keV(x_me, p = 1) = x_me * (ELECTRONMASS * 1.0e-3)^p
_me2MeV(x_me, p = 1) = x_me * (ELECTRONMASS * 1.0e-6)^p
