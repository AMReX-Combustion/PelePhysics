"""Utility functions used across ceptr."""

import ceptr.constants as cc


def qss_sorted_phase_space(mechanism, species_info, reagents):
    phi = []
    dict_species = {v: i for i, v in enumerate(species_info.all_species_list)}
    sorted_reagents = sorted(reagents.keys(), key=lambda v: dict_species[v])
    n_species = species_info.n_species
    for symbol in sorted_reagents:
        coefficient = reagents[symbol]
        if symbol in species_info.qss_species_list:
            if float(coefficient) == 1.0:
                conc = "sc_qss[%d]" % (
                    species_info.ordered_idx_map[symbol] - n_species
                )
            else:
                conc = "pow(sc_qss[%d], %f)" % (
                    species_info.ordered_idx_map[symbol] - n_species,
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


def phase_space_units(reagents):
    dim = 0.0
    for _, coefficient in reagents.items():
        dim += float(coefficient)
    return dim


def prefactor_units(units, exponent):
    return units**exponent / cc.ureg.second


def activation_energy_units():
    return cc.ureg.cal / cc.ureg.mole


def is_remove_forward(reaction_info, idx):
    return idx in reaction_info.reacRemoveIDList


def fkc_conv_inv(self, mechanism, reaction):
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


def kc_conv(mechanism, reaction):
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


def sorted_kc(mechanism, species_info, reaction):
    conv = kc_conv(mechanism, reaction)
    exparg = sorted_kc_exp_arg(mechanism, species_info, reaction)
    if conv:
        return conv + " * exp(" + exparg + ")"
    else:
        return "exp(" + exparg + ")"


def sorted_kc_exp_arg(mechanism, species_info, reaction):
    terms = []
    for _ in range(species_info.n_species):
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
            i = species_info.ordered_idx_map[symbol] - species_info.n_species
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
            i = species_info.ordered_idx_map[symbol] - species_info.n_species
            terms_qss[i] += "%sg_RT_qss[%d]" % (factor, i)
        else:
            i = species_info.ordered_idx_map[symbol]
            terms[i] += "%sg_RT[%d]" % (factor, i)

    dG = ""
    for i in range(species_info.n_species):
        if terms[i]:
            dG += terms[i]
    for i in range(species_info.nQSSspecies):
        if terms_qss[i]:
            dG += terms_qss[i]

    if dG[0:3] == " + ":
        return dG[3:]
    else:
        return "-" + dG[3:]
