#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from __future__ import absolute_import, division

import pyre
from pyre.handbook.units.energy import J, cal, kcal, kJ
from pyre.handbook.units.length import cm
from pyre.handbook.units.SI import kelvin, meter, mole, second


def equilibrium(mixture, reaction):

    dim = 0
    reagents = []

    for species, coefficient in reaction.reactants:
        reagents.append((mixture.find(species), -coefficient))
        dim -= coefficient

    for species, coefficient in reaction.products:
        reagents.append((mixture.find(species), coefficient))
        dim += coefficient

    from .Equilibrium import Equilibrium

    return Equilibrium(reagents, dim)


def forwardRate(mixture, reaction):

    # adjust the unist of the activation energy
    # the prefactor adjustment depends on the phase space dimensions and the third body
    # this logic can be found in the sections for each reaction type

    arrhenius = reaction.arrhenius
    energyUC = activationEnergyUnits(reaction.units["activation"])
    ps = phaseSpaceFactor(mixture, reaction.reactants)

    low = reaction.low
    thirdBody = reaction.thirdBody
    exponent = -dimPhaseSpace(reaction.reactants)

    if not thirdBody or low:
        exponent += 1

    prefactorUC = prefactorUnitConversion(
        reaction.units["prefactor"], exponent
    )

    # perform the unit conversion
    A, beta, E = arrhenius
    E *= energyUC
    A *= prefactorUC / kelvin ** beta
    arrhenius = (A, beta, E)

    if not low:
        lt = reaction.lt
        if lt:
            rate = pyre.models.kinetics.landau_teller(arrhenius, lt)
        else:
            rate = pyre.models.kinetics.arrhenius(arrhenius)

        if not thirdBody:
            from .Arrhenius import Arrhenius

            return Arrhenius(ps, rate)

        phaseEnhancement = enhancement(mixture, reaction)
        from .ThirdBody import ThirdBody

        return ThirdBody(ps, phaseEnhancement, rate)

    # the rest are all falloff reactions
    # adjust the units of the low coefficients
    exponent -= 1  # Pr = [M] A_o/Ai_nfty
    prefactorUC = prefactorUnitConversion(
        reaction.units["prefactor"], exponent
    )
    A, beta, E = low
    E *= energyUC
    A *= prefactorUC / kelvin ** beta
    low = (A, beta, E)
    phaseEnhancement = enhancement(mixture, reaction)

    sri = reaction.sri
    troe = reaction.troe
    if troe:
        troe = troe[0:1] + tuple([p * kelvin for p in troe[1:]])
        rate = pyre.models.kinetics.troe(arrhenius, low, troe)
    elif sri:
        adj = [sri[0], sri[1] * kelvin, sri[2] / kelvin]
        if len(sri) > 3:
            adj += [sri[3] / kelvin ** sri[4], sri[4]]
        rate = pyre.models.kinetics.sri(arrhenius, low, adj)
    else:
        rate = pyre.models.kinetics.lindemann(arrhenius, low)

    from .Falloff import Falloff

    return Falloff(ps, phaseEnhancement, rate)


def reverseRate(mixture, reaction, equilibriumCalculator, forwardRate):
    ps = phaseSpaceFactor(mixture, reaction.products)

    derived = not reaction.rev
    reversible = reaction.reversible

    if reversible and derived:
        from .ReverseRate import ReverseRate

        return ReverseRate(ps, forwardRate, equilibriumCalculator)

    if not reversible:
        exponent = 1 - dimPhaseSpace(reaction.products)
        prefactorUC = prefactorUnitConversion(
            reaction.units["prefactor"], exponent
        )
        rate = pyre.models.kinetics.arrhenius((0 * prefactorUC, 0, 0))

        from .Arrhenius import Arrhenius

        return Arrhenius(ps, rate)

    rev = reaction.rev
    pyre.debug.Firewall.verify(rev, "unknown reaction type")

    ps = phaseSpaceFactor(mixture, reaction.products)
    energyUC = activationEnergyUnits(reaction.units["activation"])

    thirdBody = reaction.thirdBody
    exponent = -dimPhaseSpace(reaction.products)
    if not thirdBody:
        exponent += 1

    prefactorUC = prefactorUnitConversion(
        reaction.units["prefactor"], exponent
    )

    A, beta, E = rev
    E *= energyUC
    A *= prefactorUC / kelvin ** beta
    rev = (A, beta, E)

    rlt = reaction.rlt
    if rlt:
        rate = pyre.models.kinetics.landau_teller(rev, rlt)
    else:
        rate = pyre.models.kinetics.arrhenius(reaction.rev)

    if not thirdBody:
        from .Arrhenius import Arrhenius

        return Arrhenius(ps, rate)

    phaseEnhancement = enhancement(mixture, reaction)
    from .ThirdBody import ThirdBody

    return ThirdBody(ps, phaseEnhancement, rate)


def activationEnergyUnits(code):
    if code == "cal/mole":
        units = cal / mole
    elif code == "kcal/mole":
        units = kcal / mole
    elif code == "joules/mole":
        units = J / mole
    elif code == "kjoules/mole":
        units = kJ / mole
    elif code == "kelvins":
        units = gas_constant * kelvin
    else:
        pyre.debug.Firewall.hit("unknown activation energy units '%s'" % code)
        return 1

    return units


def prefactorUnitConversion(code, exponent):

    if code == "mole/cm**3":
        units = mole / cm ** 3
    elif code == "moles":
        units = mole / cm ** 3
    elif code == "molecules":
        avogadro = pyre.handbook.constants.fundamental.avogadro
        units = 1.0 / avogadro / cm ** 3
    else:
        pyre.debug.Firewall.hit("unknown prefactor units '%s'" % code)
        return 1

    return units ** exponent / second


def dimPhaseSpace(reagents):
    exponent = 0
    for species, coefficient in reagents:
        exponent += coefficient

    return exponent


def phaseSpaceFactor(mixture, reagents):
    from .PhaseSpace import PhaseSpace

    species = []
    coefficients = []

    for symbol, coefficient in reagents:
        species.append(mixture.find(symbol))
        coefficients.append(coefficient)

    return PhaseSpace(species, coefficients)


def enhancement(mixture, reaction):

    pyre.debug.Firewall.verify(
        reaction.thirdBody,
        "enhancement called for reaction without third body",
    )

    species, coefficient = reaction.thirdBody
    efficiencies = reaction.efficiencies
    if not efficiencies:
        from .ThirdBodyEnhancement import ThirdBodyEnhancement

        return ThidBodyEnhancement(mixture.find(species))

    from .EfficiencyEnhancement import EfficiencyEnhancement

    mix = mixture.find("<mixture>")
    resolvedEfficiencies = []
    for symbol, factor in efficiencies:
        resolvedEfficiencies.append((mixture.find(symbol), factor - 1))

    return EfficiencyEnhancement(mix, resolvedEfficiencies)


# version
__id__ = "$Id$"

#
# End of file
