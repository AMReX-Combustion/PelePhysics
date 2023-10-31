"""Thermodynamics functions."""

import io
import itertools
from collections import OrderedDict

import numpy as np

import ceptr.writer as cw


def thermo(fstream, mechanism, species_info, syms=None):
    """Write thermodynamics routines."""
    models = analyze_thermodynamics(mechanism, species_info.nonqssa_species_list)
    if species_info.n_qssa_species > 0:
        qss_species_coeffs = analyze_thermodynamics_old(mechanism, species_info, 1)

    cv(fstream, species_info, models)
    cp(fstream, species_info, models)
    gibbs(fstream, species_info, models, 0, syms)
    if species_info.n_qssa_species > 0:
        gibbs_old(fstream, species_info, qss_species_coeffs, 1, syms)
    helmholtz(fstream, species_info, models)
    species_internal_energy(fstream, species_info, models)
    species_enthalpy(fstream, species_info, models, 0, syms)
    if species_info.n_qssa_species > 0:
        species_enthalpy_old(fstream, species_info, qss_species_coeffs, 1, syms)
    species_entropy(fstream, species_info, models)
    dcvpdtemp(fstream, species_info, models)


def generator_attributes(name):
    """Return the attributes for a thermo function generator."""
    dct = {
        "cv_R": {
            "function": cv_nasa7,
            "T4": True,
            "inv_temp": False,
            "inv_temp2": False,
            "log_temp": False,
        },
        "cp_R": {
            "function": cp_nasa7,
            "T4": True,
            "inv_temp": False,
            "inv_temp2": False,
            "log_temp": False,
        },
        "gibbs": {
            "function": gibbs_nasa7,
            "T4": True,
            "inv_temp": 1,
            "inv_temp2": False,
            "log_temp": True,
        },
        "gibbs_qss": {
            "function": gibbs_nasa7,
            "T4": True,
            "inv_temp": True,
            "inv_temp2": False,
            "log_temp": True,
        },
        "helmholtz": {
            "function": helmholtz_nasa7,
            "T4": True,
            "inv_temp": True,
            "inv_temp2": False,
            "log_temp": True,
        },
        "speciesInternalEnergy": {
            "function": internal_energy,
            "T4": True,
            "inv_temp": True,
            "inv_temp2": False,
            "log_temp": False,
        },
        "speciesEnthalpy": {
            "function": enthalpy_nasa7,
            "T4": True,
            "inv_temp": True,
            "inv_temp2": False,
            "log_temp": False,
        },
        "speciesEnthalpy_qss": {
            "function": enthalpy_nasa7,
            "T4": True,
            "inv_temp": True,
            "inv_temp2": False,
            "log_temp": False,
        },
        "speciesEntropy": {
            "function": entropy_nasa7,
            "T4": True,
            "inv_temp": False,
            "inv_temp2": False,
            "log_temp": True,
        },
        "dcvpRdT": {
            "function": dcpdtemp_nasa7,
            "T4": False,
            "inv_temp": False,
            "inv_temp2": False,
            "log_temp": False,
        },
    }
    return dct[name]


def analyze_thermodynamics_old(mechanism, species_info, qss_flag):
    """Extract information from the thermodynamics model."""
    midpoints = OrderedDict()

    spec_list = (
        species_info.qssa_species_list
        if qss_flag
        else species_info.nonqssa_species_list
    )
    for symbol in spec_list:
        species = mechanism.species(symbol)
        model = species.thermo

        assert model.n_coeffs == 15, "Unsupported thermo model."

        mid = model.coeffs[0]
        high_range = model.coeffs[1:8]
        low_range = model.coeffs[8:15]

        midpoints.setdefault(mid, []).append((species, low_range, high_range))

    return midpoints


def analyze_thermodynamics(mechanism, species_list):
    """Extract information from the thermodynamics model."""
    models = []
    for symbol in species_list:
        species = mechanism.species(symbol)
        model = species.thermo
        dct = {"species": species}

        # for nasa7
        assert model.n_coeffs == 15, "Unsupported thermo model."
        interval = []
        coeffs = [model.coeffs[1:8]]
        if not np.allclose(model.coeffs[1:8], model.coeffs[8:15], atol=1e-22):
            interval = [model.coeffs[0]]
            coeffs = [model.coeffs[8:15], model.coeffs[1:8]]

        dct["interval"] = interval
        dct["coefficients"] = coeffs
        models.append(dct)

    return models


def generate_thermo_routine(
    fstream,
    species_info,
    name,
    models,
    qss_flag,
    syms=None,
    inline=False,
):
    """Write a thermodynamics routine."""
    gen = generator_attributes(name)
    expression_generator = gen["function"]

    if not inline:
        cw.writer(
            fstream,
            f"AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void {name}(amrex::Real"
            " * species, const amrex::Real T)",
        )

    # if syms:
    #     syms_g_rt = name == "gibbs"
    #     syms_g_rt_qss = name == "gibbs_qss"
    #     syms_h_rt = name == "speciesEnthalpy"
    #     syms_h_rt_qss = name == "speciesEnthalpy_qss"

    if not inline:
        cw.writer(fstream, "{")
    cw.writer(fstream, "const amrex::Real T2 = T*T;")
    cw.writer(fstream, "const amrex::Real T3 = T*T*T;")
    if gen["T4"]:
        cw.writer(fstream, "const amrex::Real T4 = T*T*T*T;")
    if gen["inv_temp"]:
        cw.writer(fstream, "const amrex::Real invT = 1.0 / T;")
    if gen["inv_temp2"]:
        cw.writer(fstream, "const amrex::Real invT2 = invT*invT;")
    if gen["log_temp"]:
        cw.writer(fstream, "const amrex::Real logT = log(T);")
    cw.writer(fstream)

    intervals = sorted([x["interval"] for x in models])
    intervals = list(intervals for intervals, _ in itertools.groupby(intervals))

    for interval in intervals:
        cw.writer(fstream)
        if len(interval) == 0:
            cw.writer(
                fstream,
                cw.comment("species with no change across T"),
            )
        elif len(interval) == 1:
            cw.writer(
                fstream,
                cw.comment(f"species with midpoint at T={interval[0]:g} kelvin"),
            )
        for k in range(len(interval) + 1):
            if len(interval) == 1:
                if k == 0:
                    cw.writer(
                        fstream,
                        f"""if (T < {interval[0]:g}) {{""",
                    )
                else:
                    cw.writer(
                        fstream,
                        "else {",
                    )

            for model in [x for x in models if x["interval"] == interval]:
                species = model["species"]

                index = (
                    species_info.ordered_idx_map[species.name] - species_info.n_species
                    if qss_flag
                    else species_info.ordered_idx_map[species.name]
                )
                cw.writer(fstream, cw.comment(f"species {index}: {species.name}"))
                cw.writer(
                    fstream,
                    (f"result += y[{index}] * (" if inline else f"species[{index}] ="),
                )
                expression_generator(fstream, model["coefficients"][k])
                if inline:
                    spec_idx = species_info.ordered_idx_map[species.name]
                    sp = species_info.nonqssa_species[spec_idx]
                    imw = 1.0 / sp.weight
                    cw.writer(fstream, f")* {imw:.16f}")
                cw.writer(fstream, ";")

            if len(interval) == 1:
                cw.writer(
                    fstream,
                    "}",
                )

    if not inline:
        cw.writer(fstream, "}")


def generate_thermo_routine_old(
    fstream,
    species_info,
    name,
    expression_generator,
    species_coeffs,
    qss_flag,
    needs_inv_temp=0,
    needs_log_temp=False,
    syms=None,
    inline=False,
):
    """Write a thermodynamics routine."""
    midpoints = species_coeffs

    if not inline:
        cw.writer(
            fstream,
            f"AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void {name}(amrex::Real"
            " * species, const amrex::Real T)",
        )

    syms_g_rt = False
    syms_g_rt_qss = False
    syms_h_rt = False
    syms_h_rt_qss = False
    if name == "gibbs" and not (syms is None):
        syms_g_rt = True
    if name == "gibbs_qss" and not (syms is None):
        syms_g_rt_qss = True
    if name == "speciesEnthalpy" and not (syms is None):
        syms_h_rt = True
    if name == "speciesEnthalpy_qss" and not (syms is None):
        syms_h_rt_qss = True

    if not inline:
        cw.writer(fstream, "{")

    cw.writer(fstream, "const amrex::Real T2 = T*T;")
    cw.writer(fstream, "const amrex::Real T3 = T*T*T;")
    cw.writer(fstream, "const amrex::Real T4 = T*T*T*T;")
    if needs_inv_temp != 0:
        cw.writer(fstream, "const amrex::Real invT = 1.0 / T;")
    if needs_inv_temp == 2:
        cw.writer(fstream, "const amrex::Real invT2 = invT*invT;")
    if needs_log_temp:
        cw.writer(fstream, "const amrex::Real logT = log(T);")

    for mid_temp, species_list in list(midpoints.items()):
        lostream = io.StringIO()
        for species, low_range, _ in species_list:
            index = (
                species_info.ordered_idx_map[species.name] - species_info.n_species
                if qss_flag
                else species_info.ordered_idx_map[species.name]
            )
            cw.writer(lostream, cw.comment(f"species {index}: {species.name}"))
            cw.writer(
                lostream,
                (f"result += y[{index}] * (" if inline else f"species[{index}] ="),
            )
            if syms_g_rt:
                syms.g_RT_smp_tmp[mid_temp]["m"][index] = expression_generator(
                    lostream, low_range, syms
                )
            elif syms_g_rt_qss:
                syms.g_RT_qss_smp_tmp[mid_temp]["m"][index] = expression_generator(
                    lostream, low_range, syms
                )
            elif syms_h_rt:
                syms.h_RT_smp_tmp[mid_temp]["m"][index] = expression_generator(
                    lostream, low_range, syms
                )
            elif syms_h_rt_qss:
                syms.h_RT_qss_smp_tmp[mid_temp]["m"][index] = expression_generator(
                    lostream, low_range, syms
                )
            else:
                expression_generator(lostream, low_range)
            if inline:
                spec_idx = species_info.ordered_idx_map[species.name]
                sp = species_info.nonqssa_species[spec_idx]
                imw = 1.0 / sp.weight
                cw.writer(lostream, f")* {imw:.16f}")
            cw.writer(lostream, ";")

        histream = io.StringIO()
        for species, _, high_range in species_list:
            index = (
                species_info.ordered_idx_map[species.name] - species_info.n_species
                if qss_flag
                else species_info.ordered_idx_map[species.name]
            )
            cw.writer(
                histream,
                cw.comment(f"species {index}: {species.name}"),
            )
            cw.writer(
                histream,
                (f"result += y[{index}] * (" if inline else f"species[{index}] ="),
            )
            if syms_g_rt:
                syms.g_RT_smp_tmp[mid_temp]["p"][index] = expression_generator(
                    histream, high_range, syms
                )
            elif syms_g_rt_qss:
                syms.g_RT_qss_smp_tmp[mid_temp]["p"][index] = expression_generator(
                    histream, high_range, syms
                )
            elif syms_h_rt:
                syms.h_RT_smp_tmp[mid_temp]["p"][index] = expression_generator(
                    histream, high_range, syms
                )
            elif syms_h_rt_qss:
                syms.h_RT_qss_smp_tmp[mid_temp]["p"][index] = expression_generator(
                    histream, high_range, syms
                )
            else:
                expression_generator(histream, high_range)
            if inline:
                spec_idx = species_info.ordered_idx_map[species.name]
                sp = species_info.nonqssa_species[spec_idx]
                imw = 1.0 / sp.weight
                cw.writer(histream, f")* {imw:.16f}")
            cw.writer(histream, ";")

        lostr = lostream.getvalue().rstrip("\n")
        histr = histream.getvalue().rstrip("\n")
        cw.writer(fstream, "")
        if histr == lostr:
            cw.writer(
                fstream,
                cw.comment("species with no change at a midpoint T"),
            )
            cw.writer(fstream, lostr)
        else:
            cw.writer(
                fstream,
                cw.comment(f"species with midpoint at T={mid_temp:g} kelvin"),
            )
            cw.writer(
                fstream,
                f"""if (T < {mid_temp:g}) {{\n{lostr}}} else {{\n{histr}}}""",
            )
        lostream.close()
        histream.close()

    if not inline:
        cw.writer(fstream, "}")


def cv(fstream, species_info, models):
    """Write cv."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute Cv/R at the given temperature"))
    generate_thermo_routine(fstream, species_info, "cv_R", models, 0)


def cp(fstream, species_info, models):
    """Write cp."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute Cp/R at the given temperature"))
    generate_thermo_routine(fstream, species_info, "cp_R", models, 0)


def gibbs(fstream, species_info, models, qss_flag, syms=None):
    """Write Gibbs."""
    name = "gibbs_qss" if qss_flag else "gibbs"
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the g/(RT) at the given temperature"))
    generate_thermo_routine(fstream, species_info, name, models, qss_flag, syms)


def gibbs_old(fstream, species_info, species_coeffs, qss_flag, syms=None):
    """Write Gibbs."""
    name = "gibbs_qss" if qss_flag else "gibbs"
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the g/(RT) at the given temperature"))
    generate_thermo_routine_old(
        fstream,
        species_info,
        name,
        gibbs_nasa7,
        species_coeffs,
        qss_flag,
        1,
        True,
        syms,
    )


def helmholtz(fstream, species_info, models):
    """Write Helmholtz."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the a/(RT) at the given temperature"))
    generate_thermo_routine(
        fstream,
        species_info,
        "helmholtz",
        models,
        0,
    )


def species_internal_energy(fstream, species_info, models):
    """Write species internal energy."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the e/(RT) at the given temperature"))
    generate_thermo_routine(
        fstream,
        species_info,
        "speciesInternalEnergy",
        models,
        0,
    )


def species_enthalpy(fstream, species_info, models, qss_flag, syms=None):
    """Write species enthalpy."""
    name = "speciesEnthalpy_qss" if qss_flag else "speciesEnthalpy"
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("compute the h/(RT) at the given temperature (Eq 20)"),
    )

    generate_thermo_routine(
        fstream,
        species_info,
        name,
        models,
        qss_flag,
        syms,
    )


def species_enthalpy_old(fstream, species_info, species_coeffs, qss_flag, syms=None):
    """Write species enthalpy."""
    name = "speciesEnthalpy_qss" if qss_flag else "speciesEnthalpy"
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("compute the h/(RT) at the given temperature (Eq 20)"),
    )

    generate_thermo_routine_old(
        fstream,
        species_info,
        name,
        enthalpy_nasa7,
        species_coeffs,
        qss_flag,
        1,
        False,
        syms,
    )


def species_entropy(fstream, species_info, models):
    """Write species entropy."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the S/R at the given temperature (Eq 21)"))
    generate_thermo_routine(
        fstream,
        species_info,
        "speciesEntropy",
        models,
        0,
    )


def dcvpdtemp(fstream, species_info, models):
    """Write gradient of cp/cv wrt temperature."""
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature"),
    )
    generate_thermo_routine(fstream, species_info, "dcvpRdT", models, 0)


def dcpdtemp_nasa7(fstream, parameters):
    """Write NASA7 polynomial for dcpdtemp."""
    cw.writer(fstream, f"{parameters[1]:+15.8e}")
    cw.writer(fstream, f"{(parameters[2] * 2.0):+15.8e} * T")
    cw.writer(fstream, f"{(parameters[3] * 3.0):+15.8e} * T2")
    cw.writer(fstream, f"{(parameters[4] * 4.0):+15.8e} * T3")


def cv_nasa7(fstream, parameters):
    """Write NASA7 polynomial for cv."""
    cw.writer(fstream, f"{(parameters[0] - 1.0):+15.8e}")
    cw.writer(fstream, f"{parameters[1]:+15.8e} * T")
    cw.writer(fstream, f"{parameters[2]:+15.8e} * T2")
    cw.writer(fstream, f"{parameters[3]:+15.8e} * T3")
    cw.writer(fstream, f"{parameters[4]:+15.8e} * T4")


def eval_cv_species(mechanism, species, temp):
    """Evaluate cv."""
    model = mechanism.species(species.name).thermo
    return eval_cv_species_nasa7(model, temp)


def eval_cv_species_nasa7(model, temp):
    """Evaluate cv with NASA7 polynomial."""
    assert model.n_coeffs == 15, "Unsupported thermo model."
    mid = model.coeffs[0]
    high_range = model.coeffs[1:8]
    low_range = model.coeffs[8:15]

    if temp < mid:
        parameters = low_range
    else:
        parameters = high_range

    return (
        (parameters[0] - 1.0)
        + parameters[1] * temp
        + parameters[2] * temp * temp
        + parameters[3] * temp * temp * temp
        + parameters[4] * temp * temp * temp * temp
    )


def cp_nasa7(fstream, parameters):
    """Write NASA7 polynomial for cp."""
    cw.writer(fstream, f"{parameters[0]:+15.8e}")
    cw.writer(fstream, f"{parameters[1]:+15.8e} * T")
    cw.writer(fstream, f"{parameters[2]:+15.8e} * T2")
    cw.writer(fstream, f"{parameters[3]:+15.8e} * T3")
    cw.writer(fstream, f"{parameters[4]:+15.8e} * T4")


def gibbs_nasa7(fstream, parameters, syms=None):
    """Write NASA7 polynomial for Gibbs."""
    cw.writer(fstream, f"{parameters[5]:+20.15e} * invT")
    cw.writer(fstream, f"{(parameters[0] - parameters[6]):+20.15e}")
    cw.writer(fstream, f"{(-parameters[0]):+20.15e} * logT")
    cw.writer(fstream, f"{((-parameters[1] / 2)):+20.15e} * T")
    cw.writer(fstream, f"{((-parameters[2] / 6)):+20.15e} * T2")
    cw.writer(fstream, f"{((-parameters[3] / 12)):+20.15e} * T3")
    cw.writer(fstream, f"{((-parameters[4] / 20)):+20.15e} * T4")

    if syms:
        return (
            parameters[5] * syms.invT_smp
            + parameters[0]
            - parameters[6]
            + (-parameters[0]) * syms.logT_smp
            + (-parameters[1] / 2) * syms.T_smp
            + (-parameters[2] / 6) * syms.T2_smp
            + (-parameters[3] / 12) * syms.T3_smp
            + (-parameters[4] / 20) * syms.T4_smp
        )


def helmholtz_nasa7(fstream, parameters, syms=None):
    """Write NASA7 polynomial for Helmholtz."""
    cw.writer(fstream, f"{parameters[5]:+15.8e} * invT")
    cw.writer(fstream, f"{(parameters[0] - parameters[6] - 1.0):+15.8e}")
    cw.writer(fstream, f"{(-parameters[0]):+15.8e} * logT")
    cw.writer(fstream, f"{((-parameters[1] / 2)):+15.8e} * T")
    cw.writer(fstream, f"{((-parameters[2] / 6)):+15.8e} * T2")
    cw.writer(fstream, f"{((-parameters[3] / 12)):+15.8e} * T3")
    cw.writer(fstream, f"{((-parameters[4] / 20)):+15.8e} * T4")

    if syms:
        return (
            parameters[5] * syms.invT_smp
            + parameters[0]
            - parameters[6]
            - 1.0
            + (-parameters[0]) * syms.logT_smp
            + (-parameters[1] / 2) * syms.T_smp
            + (-parameters[2] / 6) * syms.T2_smp
            + (-parameters[3] / 12) * syms.T3_smp
            + (-parameters[4] / 20) * syms.T4_smp
        )


def internal_energy(fstream, parameters, syms=None):
    """Write NASA7 polynomial for internal energy."""
    cw.writer(fstream, f"{(parameters[0] - 1.0):+15.8e}")
    cw.writer(fstream, f"{((parameters[1] / 2)):+15.8e} * T")
    cw.writer(fstream, f"{((parameters[2] / 3)):+15.8e} * T2")
    cw.writer(fstream, f"{((parameters[3] / 4)):+15.8e} * T3")
    cw.writer(fstream, f"{((parameters[4] / 5)):+15.8e} * T4")
    cw.writer(fstream, f"{(parameters[5]):+15.8e} * invT")
    if syms:
        return (
            parameters[0]
            - 1.0
            + (parameters[1] / 2) * syms.T_smp
            + (parameters[2] / 3) * syms.T2_smp
            + (parameters[3] / 4) * syms.T3_smp
            + (parameters[4] / 5) * syms.T4_smp
            + (parameters[5]) * syms.invT_smp
        )


def enthalpy_nasa7(fstream, parameters, syms=None):
    """Write NASA7 polynomial for enthalpy."""
    cw.writer(fstream, f"{parameters[0]:+15.8e}")
    cw.writer(fstream, f"{((parameters[1] / 2)):+15.8e} * T")
    cw.writer(fstream, f"{((parameters[2] / 3)):+15.8e} * T2")
    cw.writer(fstream, f"{((parameters[3] / 4)):+15.8e} * T3")
    cw.writer(fstream, f"{((parameters[4] / 5)):+15.8e} * T4")
    cw.writer(fstream, f"{(parameters[5]):+15.8e} * invT")

    if syms:
        return (
            parameters[0]
            + (parameters[1] / 2) * syms.T_smp
            + (parameters[2] / 3) * syms.T2_smp
            + (parameters[3] / 4) * syms.T3_smp
            + (parameters[4] / 5) * syms.T4_smp
            + (parameters[5]) * syms.invT_smp
        )


def entropy_nasa7(fstream, parameters, syms=None):
    """Write NASA7 polynomial for entropy."""
    cw.writer(fstream, f"{parameters[0]:+15.8e} * logT")
    cw.writer(fstream, f"{(parameters[1]):+15.8e} * T")
    cw.writer(fstream, f"{((parameters[2] / 2)):+15.8e} * T2")
    cw.writer(fstream, f"{((parameters[3] / 3)):+15.8e} * T3")
    cw.writer(fstream, f"{((parameters[4] / 4)):+15.8e} * T4")
    cw.writer(fstream, f"{(parameters[6]):+15.8e}")

    if syms:
        return (
            (parameters[0]) * syms.logT_smp
            + (parameters[1]) * syms.T_smp
            + (parameters[2] / 2) * syms.T2_smp
            + (parameters[3] / 3) * syms.T3_smp
            + (parameters[4] / 4) * syms.T4_smp
            + parameters[6]
        )
