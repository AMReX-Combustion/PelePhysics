"""Utility functions used across ceptr."""
import copy
import sys
from collections import Counter
from math import isclose

import ceptr.constants as cc


def intersection(lst1, lst2):
    """Return intersection of two lists."""
    return list(set(lst1).intersection(lst2))


def qss_sorted_phase_space(mechanism, species_info, reaction, reagents, syms=None):
    """Get string of phase space."""
    record_symbolic_operations = True
    if syms is None:
        record_symbolic_operations = False
    if reaction.third_body:
        if len(reaction.third_body.efficiencies) == 1:
            if isclose(reaction.third_body.default_efficiency, 0.0):
                reagents = copy.deepcopy(
                    dict(
                        sum(
                            (
                                Counter(x)
                                for x in [
                                    reagents,
                                    reaction.third_body.efficiencies,
                                ]
                            ),
                            Counter(),
                        )
                    )
                )
    phi = []
    if record_symbolic_operations:
        phi_smp = []
    dict_species = {v: i for i, v in enumerate(species_info.all_species_list)}
    sorted_reagents = sorted(reagents.keys(), key=lambda v: dict_species[v])
    n_species = species_info.n_species
    for symbol in sorted_reagents:
        coefficient = reagents[symbol]
        if symbol in species_info.qssa_species_list:
            if float(coefficient) == 1.0:
                conc = f"sc_qss[{species_info.ordered_idx_map[symbol] - n_species}]"
                if record_symbolic_operations:
                    conc_smp = syms.sc_qss_smp[
                        species_info.ordered_idx_map[symbol] - n_species
                    ]
            else:
                if coefficient.is_integer():
                    conc = "*".join(
                        [f"sc_qss[{species_info.ordered_idx_map[symbol] - n_species}]"]
                        * int(coefficient)
                    )
                else:
                    conc = (
                        f"pow(sc_qss[{species_info.ordered_idx_map[symbol] - n_species}],"
                        f" {float(coefficient):f})"
                    )
                if record_symbolic_operations:
                    conc_smp = pow(
                        syms.sc_qss_smp[
                            species_info.ordered_idx_map[symbol] - n_species
                        ],
                        float(coefficient),
                    )
            phi += [conc]
            if record_symbolic_operations:
                phi_smp += [conc_smp]
        else:
            if float(coefficient) == 1.0:
                conc = f"sc[{species_info.ordered_idx_map[symbol]}]"
                if record_symbolic_operations:
                    conc_smp = syms.sc_smp[species_info.ordered_idx_map[symbol]]
            else:
                if float(coefficient) == 2.0:
                    conc = (
                        f"(sc[{species_info.ordered_idx_map[symbol]}] *"
                        f" sc[{species_info.ordered_idx_map[symbol]}])"
                    )
                    if record_symbolic_operations:
                        conc_smp = (
                            syms.sc_smp[species_info.ordered_idx_map[symbol]]
                            * syms.sc_smp[species_info.ordered_idx_map[symbol]]
                        )
                else:
                    conc = (
                        f"pow(sc[{species_info.ordered_idx_map[symbol]}],"
                        f" {float(coefficient):f})"
                    )
                    if record_symbolic_operations:
                        conc_smp = pow(
                            syms.sc_smp[species_info.ordered_idx_map[symbol]],
                            float(coefficient),
                        )
            phi += [conc]
            if record_symbolic_operations:
                phi_smp += [conc_smp]
    if not record_symbolic_operations:
        return "*".join(phi)
    else:
        if syms.remove_1:
            out_smp = 1
        else:
            out_smp = 1.0
        for phi_val_smp in phi_smp:
            out_smp *= phi_val_smp
        out_smp = syms.convert_symb_to_int(out_smp)
        return "*".join(phi), out_smp


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
    record_symbolic_operations = True
    if syms is None:
        record_symbolic_operations = False
    dim = 0
    if record_symbolic_operations:
        dim_smp = 0
    for _, coefficient in reaction.reactants.items():
        dim -= coefficient
        if record_symbolic_operations:
            dim_smp -= coefficient
    # flip the signs
    for _, coefficient in reaction.products.items():
        dim += coefficient
        if record_symbolic_operations:
            dim_smp += coefficient

    conversion_smp = 1.0
    if record_symbolic_operations and syms.remove_1:
        conversion_smp = 1
    if dim == 0:
        conversion = ""
    elif dim > 0:
        if dim == 1.0:
            conversion = "*".join(["refCinv"])
            if record_symbolic_operations:
                conversion_smp *= syms.refCinv_smp
        else:
            if dim.is_integer():
                mult_str = "*".join(["refCinv"] * int(dim))
                conversion = "*".join([f"({mult_str})"])
                if record_symbolic_operations:
                    conversion_smp *= syms.refCinv_smp * syms.refCinv_smp
            else:
                conversion = "*".join([f"pow(refCinv, {dim:f})"])
                if record_symbolic_operations:
                    conversion_smp *= syms.refCinv_smp**dim
    else:
        if dim == -1.0:
            conversion = "*".join(["refC"])
            if record_symbolic_operations:
                conversion_smp *= syms.refC_smp
        else:
            conversion = "*".join([f"pow(refC, {abs(dim):f})"])
            if record_symbolic_operations:
                conversion_smp *= syms.refC_smp**dim

    if record_symbolic_operations:
        conversion_smp = syms.convert_symb_to_int(conversion_smp)
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
                conversion = "*".join([f"pow(refC,{dim:f})"])
    else:
        if dim == -1.0:
            conversion = "*".join(["refCinv"])
        else:
            conversion = "*".join([f"pow(refCinv,{abs(dim):f})"])

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
    record_symbolic_operations = True
    if syms is None:
        record_symbolic_operations = False
    terms = []
    if record_symbolic_operations:
        terms_smp = []
    for _ in range(species_info.n_species):
        terms.append("")
        if record_symbolic_operations:
            terms_smp.append(0.0)
    terms_qss = []
    if record_symbolic_operations:
        terms_qss_smp = []
    for _ in range(species_info.n_qssa_species):
        terms_qss.append("")
        if record_symbolic_operations:
            terms_qss_smp.append(0.0)

    for symbol, coefficient in reaction.reactants.items():
        if coefficient == 1.0:
            factor = " + "
            if record_symbolic_operations:
                factor_smp = 1.0
        else:
            factor = f" + {coefficient:f}*"
            if record_symbolic_operations:
                factor_smp = coefficient

        if symbol in species_info.qssa_species_list:
            i = species_info.ordered_idx_map[symbol] - species_info.n_species
            terms_qss[i] += f"{factor}g_RT_qss[{i}]"
            if record_symbolic_operations:
                factor_smp = syms.convert_symb_to_int(factor_smp)
                terms_qss_smp[i] += factor_smp * syms.g_RT_qss_smp[i]
        else:
            i = species_info.ordered_idx_map[symbol]
            terms[i] += f"{factor}g_RT[{i}]"
            if record_symbolic_operations:
                factor_smp = syms.convert_symb_to_int(factor_smp)
                terms_smp[i] += factor_smp * syms.g_RT_smp[i]

    for symbol, coefficient in reaction.products.items():
        if coefficient == 1.0:
            factor = " - "  # flip the signs
            if record_symbolic_operations:
                factor_smp = -1.0
        else:
            factor = f" - {coefficient:f}*"
            if record_symbolic_operations:
                factor_smp = -coefficient

        if symbol in species_info.qssa_species_list:
            i = species_info.ordered_idx_map[symbol] - species_info.n_species
            terms_qss[i] += f"{factor}g_RT_qss[{i}]"
            if record_symbolic_operations:
                factor_smp = syms.convert_symb_to_int(factor_smp)
                terms_qss_smp[i] += factor_smp * syms.g_RT_qss_smp[i]
        else:
            i = species_info.ordered_idx_map[symbol]
            terms[i] += f"{factor}g_RT[{i}]"
            if record_symbolic_operations:
                factor_smp = syms.convert_symb_to_int(factor_smp)
                terms_smp[i] += factor_smp * syms.g_RT_smp[i]

    dg = ""
    if record_symbolic_operations:
        dg_smp = 0
    for i in range(species_info.n_species):
        if terms[i]:
            dg += terms[i]
            if record_symbolic_operations:
                dg_smp += terms_smp[i]
    for i in range(species_info.n_qssa_species):
        if terms_qss[i]:
            dg += terms_qss[i]
            if record_symbolic_operations:
                dg_smp += terms_qss_smp[i]

    if dg[0:3] == " + ":
        # print("p dg = ", dg)
        # print("p dg_smp = ", dg_smp)
        if record_symbolic_operations:
            dg_smp = syms.convert_symb_to_int(dg_smp)
            return dg[3:], dg_smp
        else:
            return dg[3:]
    else:
        # print("m dg = ", dg)
        # print("m dg_smp = ", dg_smp)
        if record_symbolic_operations:
            dg_smp = syms.convert_symb_to_int(dg_smp)
            return "-" + dg[3:], dg_smp
        else:
            return "-" + dg[3:]


def enhancement_d(mechanism, species_info, reaction, syms=None):
    """Write get enhancement."""
    record_symbolic_operations = True
    if syms is None:
        record_symbolic_operations = False

    third_body = reaction.third_body is not None
    falloff = reaction.rate.type == "falloff"
    if not third_body and not falloff:
        print("enhancement_d called for a reaction without a third body")
        sys.exit(1)

    if not reaction.third_body:
        print("FIXME EFFICIENCIES")
        sys.exit(1)
        species, coefficient = third_body
        if species == "<mixture>":
            if record_symbolic_operations:
                return "mixture", syms.mixture_smp
            else:
                return "mixture"
        if record_symbolic_operations:
            return (
                f"sc[{species_info.ordered_idx_map[species]}]",
                syms.sc_smp[species_info.ordered_idx_map[species]],
            )
        else:
            return f"sc[{species_info.ordered_idx_map[species]}]"

    efficiencies = reaction.third_body.efficiencies
    alpha = ["mixture"]
    if record_symbolic_operations:
        alpha_smp = [syms.mixture_smp]
    dict_species = {v: i for i, v in enumerate(species_info.all_species_list)}
    sorted_efficiencies = sorted(efficiencies.keys(), key=lambda v: dict_species[v])
    for symbol in sorted_efficiencies:
        efficiency = efficiencies[symbol]
        if symbol not in species_info.qssa_species_list:
            factor = f"( {efficiency:.15g} - 1)"
            if record_symbolic_operations:
                factor_smp = efficiency - 1
            if (efficiency - 1) != 0:
                conc = f"sc[{species_info.ordered_idx_map[symbol]}]"
                if record_symbolic_operations:
                    conc_smp = syms.sc_smp[species_info.ordered_idx_map[symbol]]
                if (efficiency - 1) == 1:
                    alpha.append(f"{conc}")
                    if record_symbolic_operations:
                        alpha_smp.append(conc_smp)
                else:
                    alpha.append(f"{factor}*{conc}")
                    if record_symbolic_operations:
                        factor_smp = syms.convert_symb_to_int(factor_smp)
                        alpha_smp.append(factor_smp * conc_smp)

    if record_symbolic_operations:
        enhancement_smp = 0.0
        for alpha_val in alpha_smp:
            enhancement_smp += alpha_val

    if record_symbolic_operations:
        enhancement_smp = syms.convert_symb_to_int(enhancement_smp)
        return " + ".join(alpha).replace("+ -", "- "), enhancement_smp
    else:
        return " + ".join(alpha).replace("+ -", "- ")
