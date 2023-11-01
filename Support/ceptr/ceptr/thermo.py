"""Thermodynamics functions."""

import itertools

import numpy as np

import ceptr.writer as cw


def thermo(fstream, mechanism, species_info, syms=None):
    """Write thermodynamics routines."""
    models = analyze_thermodynamics(mechanism, species_info.nonqssa_species_list)
    if species_info.n_qssa_species > 0:
        qss_models = analyze_thermodynamics(mechanism, species_info.qssa_species_list)

    cv(fstream, species_info, models)
    cp(fstream, species_info, models)
    gibbs(fstream, species_info, models, 0, syms)
    if species_info.n_qssa_species > 0:
        gibbs(fstream, species_info, qss_models, 1, syms)
    helmholtz(fstream, species_info, models)
    species_internal_energy(fstream, species_info, models)
    species_enthalpy(fstream, species_info, models, 0, syms)
    if species_info.n_qssa_species > 0:
        species_enthalpy(fstream, species_info, qss_models, 1, syms)
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
                if syms:
                    models_smp = (
                        syms.models_smp_tmp if not qss_flag else syms.models_qss_smp_tmp
                    )
                    model_smp = [
                        x for x in models_smp if x["species"].name == species.name
                    ][0]
                    model_smp[name][k] = expression_generator(
                        fstream, model["coefficients"][k], syms
                    )
                else:
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


def param2str(param, sfx="", fmt="+15.8e"):
    return f"{param:{fmt}} {sfx}" if param != 0.0 else ""


def dcpdtemp_nasa7(fstream, parameters):
    """Write NASA7 polynomial for dcpdtemp."""
    expression = (
        param2str(parameters[1])
        + param2str(parameters[2] * 2.0, "* T")
        + param2str(parameters[3] * 3.0, "* T2")
        + param2str(parameters[4] * 4.0, "* T3")
    )
    cw.writer(fstream, expression if expression else "0.0")


def cv_nasa7(fstream, parameters):
    """Write NASA7 polynomial for cv."""
    expression = (
        param2str(parameters[0] - 1.0)
        + param2str(parameters[1], "* T")
        + param2str(parameters[2], " * T2")
        + param2str(parameters[3], "* T3")
        + param2str(parameters[4], "* T4")
    )
    cw.writer(fstream, expression if expression else "0.0")


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
    expression = (
        param2str(parameters[0])
        + param2str(parameters[1], "* T")
        + param2str(parameters[2], "* T2")
        + param2str(parameters[3], "* T3")
        + param2str(parameters[4], "* T4")
    )
    cw.writer(fstream, expression if expression else "0.0")


def gibbs_nasa7(fstream, parameters, syms=None):
    """Write NASA7 polynomial for Gibbs."""
    expression = (
        param2str(parameters[5], "* invT", "+20.15e")
        + param2str(parameters[0] - parameters[6], "", "+20.15e")
        + param2str(-parameters[0], "* logT", "+20.15e")
        + param2str((-parameters[1] / 2), "* T", "+20.15e")
        + param2str((-parameters[2] / 6), "* T2", "+20.15e")
        + param2str((-parameters[3] / 12), "* T3", "+20.15e")
        + param2str((-parameters[4] / 20), "* T4", "+20.15e")
    )
    cw.writer(fstream, expression if expression else "0.0")

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
    expression = (
        param2str(parameters[5], "* invT")
        + param2str(parameters[0] - parameters[6] - 1.0, "")
        + param2str(-parameters[0], "* logT")
        + param2str((-parameters[1] / 2), "* T")
        + param2str((-parameters[2] / 6), "* T2")
        + param2str((-parameters[3] / 12), "* T3")
        + param2str((-parameters[4] / 20), "* T4")
    )
    cw.writer(fstream, expression if expression else "0.0")

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
    expression = (
        param2str(parameters[0] - 1.0, "")
        + param2str(parameters[1] / 2, "* T")
        + param2str(parameters[2] / 3, "* T2")
        + param2str(parameters[3] / 4, "* T3")
        + param2str(parameters[4] / 5, "* T4")
        + param2str(parameters[5], "* invT")
    )
    cw.writer(fstream, expression if expression else "0.0")

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
    expression = (
        param2str(parameters[0], "")
        + param2str(parameters[1] / 2, "* T")
        + param2str(parameters[2] / 3, "* T2")
        + param2str(parameters[3] / 4, "* T3")
        + param2str(parameters[4] / 5, "* T4")
        + param2str(parameters[5], "* invT")
    )
    cw.writer(fstream, expression if expression else "0.0")

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
    expression = (
        param2str(parameters[0], "* logT")
        + param2str(parameters[1], "* T")
        + param2str(parameters[2] / 2, "* T2")
        + param2str(parameters[3] / 3, "* T3")
        + param2str(parameters[4] / 4, "* T4")
        + param2str(parameters[6])
    )
    cw.writer(fstream, expression if expression else "0.0")

    if syms:
        return (
            (parameters[0]) * syms.logT_smp
            + (parameters[1]) * syms.T_smp
            + (parameters[2] / 2) * syms.T2_smp
            + (parameters[3] / 3) * syms.T3_smp
            + (parameters[4] / 4) * syms.T4_smp
            + parameters[6]
        )
