"""Utility functions used across ceptr."""
import copy
import sys
from collections import Counter

import ceptr.constants as cc


def intersection(lst1, lst2):
    """Return intersection of two lists."""
    return list(set(lst1).intersection(lst2))


def qss_sorted_phase_space(mechanism, species_info, reaction, reagents):
    """Get string of phase space."""
    if hasattr(reaction, "efficiencies"):
        if len(reaction.efficiencies) == 1:
            reagents = copy.deepcopy(
                dict(
                    sum(
                        (
                            Counter(x)
                            for x in [reagents, reaction.efficiencies]
                        ),
                        Counter(),
                    )
                )
            )

    phi = []
    dict_species = {v: i for i, v in enumerate(species_info.all_species_list)}
    sorted_reagents = sorted(reagents.keys(), key=lambda v: dict_species[v])
    n_species = species_info.n_species
    for symbol in sorted_reagents:
        coefficient = reagents[symbol]
        if symbol in species_info.qssa_species_list:
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
    """Return dimension for phase space."""
    dim = 0.0
    for _, coefficient in reagents.items():
        dim += float(coefficient)
    return dim


def prefactor_units(units, exponent):
    """Return prefactor units."""
    return units**exponent / cc.ureg.second


def activation_energy_units():
    """Return activation energy units."""
    return cc.ureg.cal / cc.ureg.mole


def is_remove_forward(reaction_info, idx):
    """Return if the reaction is in the removed list."""
    return idx in reaction_info.remove_id_list


def fkc_conv_inv(self, mechanism, reaction, syms=None):
    """Return fkc_conv_inv."""
    recordSymbolicOperations = True
    if syms is None:
        recordSymbolicOperations = False
    dim = 0
    if recordSymbolicOperations: dim_smp = 0
    for _, coefficient in reaction.reactants.items():
        dim -= coefficient
        if recordSymbolicOperations: dim_smp -= coefficient
    # flip the signs
    for _, coefficient in reaction.products.items():
        dim += coefficient
        if recordSymbolicOperations: dim_smp += coefficient

    conversion_smp = 1.0
    if dim == 0:
        conversion = ""
    elif dim > 0:
        if dim == 1.0:
            conversion = "*".join(["refCinv"])
            if recordSymbolicOperations: conversion_smp *= syms.refCinv_smp
        else:
            if dim == 2.0:
                conversion = "*".join(["(refCinv * refCinv)"])
                if recordSymbolicOperations: conversion_smp *= syms.refCinv_smp * syms.refCinv_smp
            else:
                conversion = "*".join(["pow(refCinv, %f)" % dim])
                if recordSymbolicOperations: conversion_smp *= syms.refCinv_smp**dim
    else:
        if dim == -1.0:
            conversion = "*".join(["refC"])
            if recordSymbolicOperations: conversion_smp *= syms.refC_smp
        else:
            conversion = "*".join(["pow(refC, %f)" % abs(dim)])
            if recordSymbolicOperations: conversion_smp *= syms.refC_smp ** dim

    if recordSymbolicOperations:
        return conversion, conversion_smp
    else:
        return conversion


def kc_conv(mechanism, reaction):
    """Return kc_conv."""
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
    """Return sorted kc."""
    conv = kc_conv(mechanism, reaction)
    exparg = sorted_kc_exp_arg(mechanism, species_info, reaction)
    if conv:
        return conv + " * exp(" + exparg + ")"
    else:
        return "exp(" + exparg + ")"


def sorted_kc_exp_arg(mechanism, species_info, reaction, syms=None):
    """Return sorted kc exponent argument."""
    recordSymbolicOperations = True
    if syms is None:
        recordSymbolicOperations = False
    terms = []
    if recordSymbolicOperations: terms_smp = []
    for _ in range(species_info.n_species):
        terms.append("")
        if recordSymbolicOperations: terms_smp.append(0.0)
    terms_qss = []
    if recordSymbolicOperations: terms_qss_smp = []
    for _ in range(species_info.n_qssa_species):
        terms_qss.append("")
        if recordSymbolicOperations: terms_qss_smp.append(0.0)

    for symbol, coefficient in reaction.reactants.items():
        if coefficient == 1.0:
            factor = " + "
            if recordSymbolicOperations: factor_smp = 1.0
        else:
            factor = " + %f*" % coefficient
            if recordSymbolicOperations: factor_smp = coefficient 
            
 
        if symbol in species_info.qssa_species_list:
            i = species_info.ordered_idx_map[symbol] - species_info.n_species
            terms_qss[i] += "%sg_RT_qss[%d]" % (factor, i)
            if recordSymbolicOperations: terms_qss_smp[i] +=  factor_smp * syms.g_RT_qss_smp[i]
        else:
            i = species_info.ordered_idx_map[symbol]
            terms[i] += "%sg_RT[%d]" % (factor, i)
            if recordSymbolicOperations: terms_smp[i] +=  factor_smp * syms.g_RT_smp[i]


    for symbol, coefficient in reaction.products.items():
        if coefficient == 1.0:
            factor = " - "  # flip the signs
            if recordSymbolicOperations: factor_smp = -1.0
        else:
            factor = " - %f*" % coefficient
            if recordSymbolicOperations: factor_smp = -coefficient
            

        if symbol in species_info.qssa_species_list:
            i = species_info.ordered_idx_map[symbol] - species_info.n_species
            terms_qss[i] += "%sg_RT_qss[%d]" % (factor, i)
            if recordSymbolicOperations: terms_qss_smp[i] +=  factor_smp * syms.g_RT_qss_smp[i]
        else:
            i = species_info.ordered_idx_map[symbol]
            terms[i] += "%sg_RT[%d]" % (factor, i)
            if recordSymbolicOperations: terms_smp[i] +=  factor_smp * syms.g_RT_smp[i]

    dg = ""
    if recordSymbolicOperations: dg_smp = 0
    for i in range(species_info.n_species):
        if terms[i]:
            dg += terms[i]
            if recordSymbolicOperations: dg_smp += terms_smp[i]
    for i in range(species_info.n_qssa_species):
        if terms_qss[i]:
            dg += terms_qss[i]
            if recordSymbolicOperations: dg_smp += terms_qss_smp[i]

    if dg[0:3] == " + ":
        #print("p dg = ", dg)
        #print("p dg_smp = ", dg_smp)
        if recordSymbolicOperations: 
            return dg[3:], dg_smp
        else: 
            return dg[3:]
    else:
        #print("m dg = ", dg)
        #print("m dg_smp = ", dg_smp)
        if recordSymbolicOperations: 
            return "-" + dg[3:], dg_smp
        else: 
            return "-" + dg[3:]
        

def enhancement_d(mechanism, species_info, reaction, syms=None):
    """Write get enhancement."""
    recordSymbolicOperations = True
    if syms is None:
        recordSymbolicOperations = False
   
    third_body = reaction.reaction_type == "three-body"
    falloff = reaction.reaction_type == "falloff"
    if not third_body and not falloff:
        print("enhancement_d called for a reaction without a third body")
        sys.exit(1)

    if not hasattr(reaction, "efficiencies"):
        print("FIXME EFFICIENCIES")
        sys.exit(1)
        species, coefficient = third_body
        if species == "<mixture>":
            if recordSymbolicOperations:
                return "mixture", syms.mixture_smp
            else:
                return "mixture"
        if recordSymbolicOperations:
            return "sc[%d]" % species_info.ordered_idx_map[species], syms.sc_smp[species_info.ordered_idx_map[species]]
        else:
            return "sc[%d]" % species_info.ordered_idx_map[species]

    efficiencies = reaction.efficiencies
    alpha = ["mixture"]
    if recordSymbolicOperations: alpha_smp = [syms.mixture_smp]
    dict_species = {v: i for i, v in enumerate(species_info.all_species_list)}
    sorted_efficiencies = sorted(
        efficiencies.keys(), key=lambda v: dict_species[v]
    )
    for symbol in sorted_efficiencies:
        efficiency = efficiencies[symbol]
        if symbol not in species_info.qssa_species_list:
            factor = "( %.15g - 1)" % (efficiency)
            if recordSymbolicOperations: factor_smp = efficiency - 1
            if (efficiency - 1) != 0:
                conc = "sc[%d]" % species_info.ordered_idx_map[symbol]
                if recordSymbolicOperations: conc_smp = syms.sc_smp[species_info.ordered_idx_map[symbol]]
                if (efficiency - 1) == 1:
                    alpha.append("%s" % (conc))
                    if recordSymbolicOperations: alpha_smp.append(conc_smp)
                else:
                    alpha.append("%s*%s" % (factor, conc))
                    if recordSymbolicOperations: alpha_smp.append(factor_smp*conc_smp)

    if recordSymbolicOperations: 
        enhancement_smp = 0.0
        for alphaVal in alpha_smp:
            enhancement_smp += (alphaVal)

    if recordSymbolicOperations:
        return " + ".join(alpha).replace("+ -", "- "), enhancement_smp
    else:
        return " + ".join(alpha).replace("+ -", "- ")
