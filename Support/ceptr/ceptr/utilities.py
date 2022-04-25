"""Utility functions used across ceptr."""

import ceptr.constants as cc


def QSSsortedPhaseSpace(mechanism, species_info, reagents):
    phi = []
    sorted_reagents = sorted(
        reagents,
        key=lambda x: next(
            (y for y in species_info.all_species if y.name == x[0]), None
        ).idx,
    )
    nSpecies = species_info.nSpecies
    for symbol in sorted_reagents:
        coefficient = reagents[symbol]
        if symbol in species_info.qss_species_list:
            if float(coefficient) == 1.0:
                conc = "sc_qss[%d]" % (
                    species_info.ordered_idx_map[symbol] - nSpecies
                )
            else:
                conc = "pow(sc_qss[%d], %f)" % (
                    species_info.ordered_idx_map[symbol] - nSpecies,
                    float(coefficient),
                )
            phi += [conc]
        else:
            if float(coefficient) == 1.0:
                conc = "sc[%d]" % species_info.ordered_idx_map[symbol]
            else:
                if float(coefficient) == 2.0:
                    conc = "(sc[%d] * sc[%d])" % (
                        species_info.ordered_idx_map[symbol],
                        species_info.ordered_idx_map[symbol],
                    )
                else:
                    conc = "pow(sc[%d], %f)" % (
                        species_info.ordered_idx_map[symbol],
                        float(coefficient),
                    )
            phi += [conc]
    return "*".join(phi)


def phaseSpaceUnits(reagents):
    dim = 0.0
    for _, coefficient in reagents.items():
        dim += float(coefficient)
    return dim


def prefactorUnits(units, exponent):
    return units**exponent / cc.ureg.second


def activationEnergyUnits():
    return cc.ureg.cal / cc.ureg.mole


def isRemoveForward(reaction_info, idx):
    return idx in reaction_info.reacRemoveIDList


def fKcConvInv(self, mechanism, reaction):
    dim = 0
    for _, coefficient in reaction.reactants.items():
        dim -= coefficient
    # flip the signs
    for _, coefficient in reaction.products.items():
        dim += coefficient

    if dim == 0:
        conversion = ""
    elif dim > 0:
        if dim == 1.0:
            conversion = "*".join(["refCinv"])
        else:
            if dim == 2.0:
                conversion = "*".join(["(refCinv * refCinv)"])
            else:
                conversion = "*".join(["pow(refCinv, %f)" % dim])
    else:
        if dim == -1.0:
            conversion = "*".join(["refC"])
        else:
            conversion = "*".join(["pow(refC, %f)" % abs(dim)])

    return conversion


def KcConv(mechanism, reaction):
    dim = 0
    for _, coefficient in reaction.reactants.items():
        dim -= coefficient
    # flip the signs
    for _, coefficient in reaction.products.items():
        dim += coefficient

    if dim == 0:
        conversion = ""
    elif dim > 0:
        if dim == 1.0:
            conversion = "*".join(["refC"])
        else:
            if dim == 2.0:
                conversion = "*".join(["(refC * refC)"])
            else:
                conversion = "*".join(["pow(refC,%f)" % dim])
    else:
        if dim == -1.0:
            conversion = "*".join(["refCinv"])
        else:
            conversion = "*".join(["pow(refCinv,%f)" % abs(dim)])

    return conversion


def sortedKc(mechanism, species_info, reaction):
    conv = KcConv(mechanism, reaction)
    exparg = sortedKcExpArg(mechanism, species_info, reaction)
    if conv:
        return conv + " * exp(" + exparg + ")"
    else:
        return "exp(" + exparg + ")"


def sortedKcExpArg(mechanism, species_info, reaction):
    terms = []
    for _ in range(species_info.nSpecies):
        terms.append("")
    terms_qss = []
    for _ in range(species_info.nQSSspecies):
        terms_qss.append("")

    for symbol, coefficient in reaction.reactants.items():
        if coefficient == 1.0:
            factor = " + "
        else:
            factor = " + %f*" % coefficient

        if symbol in species_info.qss_species_list:
            i = species_info.ordered_idx_map[symbol] - species_info.nSpecies
            terms_qss[i] += "%sg_RT_qss[%d]" % (factor, i)
        else:
            i = species_info.ordered_idx_map[symbol]
            terms[i] += "%sg_RT[%d]" % (factor, i)

    for symbol, coefficient in reaction.products.items():
        if coefficient == 1.0:
            factor = " - "  # flip the signs
        else:
            factor = " - %f*" % coefficient

        if symbol in species_info.qss_species_list:
            i = species_info.ordered_idx_map[symbol] - species_info.nSpecies
            terms_qss[i] += "%sg_RT_qss[%d]" % (factor, i)
        else:
            i = species_info.ordered_idx_map[symbol]
            terms[i] += "%sg_RT[%d]" % (factor, i)

    dG = ""
    for i in range(species_info.nSpecies):
        if terms[i]:
            dG += terms[i]
    for i in range(species_info.nQSSspecies):
        if terms_qss[i]:
            dG += terms_qss[i]

    if dG[0:3] == " + ":
        return dG[3:]
    else:
        return "-" + dG[3:]
