"""Thermodynamics functions."""
import sys
from collections import OrderedDict

import ceptr.writer as cw


def thermo(fstream, mechanism, species_info, syms=None):
    """Write thermodynamics routines."""
    species_coeffs = analyze_thermodynamics(mechanism, species_info, 0)
    if species_info.n_qssa_species > 0:
        qss_species_coeffs = analyze_thermodynamics(mechanism, species_info, 1)

    cv(fstream, species_info, species_coeffs)
    cp(fstream, species_info, species_coeffs)
    gibbs(fstream, species_info, species_coeffs, 0, syms)
    if species_info.n_qssa_species > 0:
        gibbs(fstream, species_info, qss_species_coeffs, 1, syms)
    helmholtz(fstream, species_info, species_coeffs)
    species_internal_energy(fstream, species_info, species_coeffs)
    species_enthalpy(fstream, species_info, species_coeffs, 0, syms)
    if species_info.n_qssa_species > 0:
        species_enthalpy(fstream, species_info, qss_species_coeffs, 1, syms)
    species_entropy(fstream, species_info, species_coeffs)


def analyze_thermodynamics(mechanism, species_info, qss_flag):
    """Extract information from the thermodynamics model."""
    low_temp = 0.0
    high_temp = 1000000.0

    midpoints = OrderedDict()

    if qss_flag:
        for symbol in species_info.qssa_species_list:
            species = mechanism.species(symbol)
            model = species.thermo

            if not model.n_coeffs == 15:
                print("Unsupported thermo model.")
                sys.exit(1)

            lo_temp = model.min_temp
            hi_temp = model.max_temp
            if low_temp < lo_temp:
                low_temp = lo_temp
            if hi_temp < high_temp:
                high_temp = hi_temp
            mid = model.coeffs[0]
            high_range = model.coeffs[1:8]
            low_range = model.coeffs[8:15]

            midpoints.setdefault(mid, []).append(
                (species, low_range, high_range)
            )

    else:
        for symbol in species_info.nonqssa_species_list:
            species = mechanism.species(symbol)
            model = species.thermo

            if not model.n_coeffs == 15:
                print("Unsupported thermo model.")
                sys.exit(1)

            lo_temp = model.min_temp
            hi_temp = model.max_temp
            if low_temp < lo_temp:
                low_temp = lo_temp
            if hi_temp < high_temp:
                high_temp = hi_temp
            mid = model.coeffs[0]
            high_range = model.coeffs[1:8]
            low_range = model.coeffs[8:15]

            midpoints.setdefault(mid, []).append(
                (species, low_range, high_range)
            )

    species_info.low_temp = low_temp
    species_info.high_temp = high_temp
    return low_temp, high_temp, midpoints


def generate_thermo_routine(
    fstream,
    species_info,
    name,
    expression_generator,
    species_coeffs,
    qss_flag,
    needs_inv_temp=0,
    syms=None,
):
    """Write a thermodynamics routine."""
    low_temp, high_temp, midpoints = species_coeffs

    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void %s(amrex::Real *"
        " species, const amrex::Real *  tc)" % name,
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

    cw.writer(fstream, "{")
    # declarations
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("temperature"))
    cw.writer(fstream, "const amrex::Real T = tc[1];")
    if needs_inv_temp != 0:
        cw.writer(fstream, "const amrex::Real invT = 1 / T;")
    if needs_inv_temp == 2:
        cw.writer(fstream, "const amrex::Real invT2 = invT*invT;")

    # if name=="gibbs_qss":
    #    print("name = ", name)
    #    for mid_temp, species_list in list(midpoints.items()):
    #        print("midpoints = ", mid_temp)
    #    stop
    # stop
    for mid_temp, species_list in list(midpoints.items()):
        cw.writer(fstream, "")
        cw.writer(
            fstream,
            cw.comment("species with midpoint at T=%g kelvin" % mid_temp),
        )
        cw.writer(fstream, "if (T < %g) {" % mid_temp)

        for species, low_range, _ in species_list:
            if qss_flag:
                cw.writer(
                    fstream,
                    cw.comment(
                        "species %d: %s"
                        % (
                            species_info.ordered_idx_map[species.name]
                            - species_info.n_species,
                            species.name,
                        )
                    ),
                )
                cw.writer(
                    fstream,
                    "species[%d] ="
                    % (
                        species_info.ordered_idx_map[species.name]
                        - species_info.n_species
                    ),
                )
            else:
                cw.writer(
                    fstream,
                    cw.comment(
                        "species %d: %s"
                        % (
                            species_info.ordered_idx_map[species.name],
                            species.name,
                        )
                    ),
                )
                cw.writer(
                    fstream,
                    "species[%d] ="
                    % (species_info.ordered_idx_map[species.name]),
                )
            if syms_g_rt:
                index = species_info.ordered_idx_map[species.name]
                syms.g_RT_smp_tmp[mid_temp]["m"][index] = expression_generator(
                    fstream, low_range, syms
                )
            elif syms_g_rt_qss:
                index = (
                    species_info.ordered_idx_map[species.name]
                    - species_info.n_species
                )
                syms.g_RT_qss_smp_tmp[mid_temp]["m"][
                    index
                ] = expression_generator(fstream, low_range, syms)
            elif syms_h_rt:
                index = species_info.ordered_idx_map[species.name]
                syms.h_RT_smp_tmp[mid_temp]["m"][index] = expression_generator(
                    fstream, low_range, syms
                )
            elif syms_h_rt_qss:
                index = (
                    species_info.ordered_idx_map[species.name]
                    - species_info.n_species
                )
                syms.h_RT_qss_smp_tmp[mid_temp]["m"][
                    index
                ] = expression_generator(fstream, low_range, syms)
            else:
                expression_generator(fstream, low_range)

        cw.writer(fstream, "} else {")

        for species, _, high_range in species_list:
            if qss_flag:
                cw.writer(
                    fstream,
                    cw.comment(
                        "species %d: %s"
                        % (
                            species_info.ordered_idx_map[species.name]
                            - species_info.n_species,
                            species.name,
                        )
                    ),
                )
                cw.writer(
                    fstream,
                    "species[%d] ="
                    % (
                        species_info.ordered_idx_map[species.name]
                        - species_info.n_species
                    ),
                )
            else:
                cw.writer(
                    fstream,
                    cw.comment(
                        "species %d: %s"
                        % (
                            species_info.ordered_idx_map[species.name],
                            species.name,
                        )
                    ),
                )
                cw.writer(
                    fstream,
                    "species[%d] ="
                    % (species_info.ordered_idx_map[species.name]),
                )
            if syms_g_rt:
                index = species_info.ordered_idx_map[species.name]
                syms.g_RT_smp_tmp[mid_temp]["p"][index] = expression_generator(
                    fstream, high_range, syms
                )
            elif syms_g_rt_qss:
                index = (
                    species_info.ordered_idx_map[species.name]
                    - species_info.n_species
                )
                syms.g_RT_qss_smp_tmp[mid_temp]["p"][
                    index
                ] = expression_generator(fstream, high_range, syms)
            elif syms_h_rt:
                index = species_info.ordered_idx_map[species.name]
                syms.h_RT_smp_tmp[mid_temp]["p"][index] = expression_generator(
                    fstream, high_range, syms
                )
            elif syms_h_rt_qss:
                index = (
                    species_info.ordered_idx_map[species.name]
                    - species_info.n_species
                )
                syms.h_RT_qss_smp_tmp[mid_temp]["p"][
                    index
                ] = expression_generator(fstream, high_range, syms)
            else:
                expression_generator(fstream, high_range)

        cw.writer(fstream, "}")

    cw.writer(fstream, "}")


def cv(fstream, species_info, species_coeffs):
    """Write cv."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute Cv/R at the given temperature"))
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    generate_thermo_routine(
        fstream, species_info, "cv_R", cv_nasa, species_coeffs, 0
    )


def cp(fstream, species_info, species_coeffs):
    """Write cp."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute Cp/R at the given temperature"))
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    generate_thermo_routine(
        fstream, species_info, "cp_R", cp_nasa, species_coeffs, 0
    )


def gibbs(fstream, species_info, species_coeffs, qss_flag, syms=None):
    """Write Gibbs."""
    if qss_flag:
        name = "gibbs_qss"
    else:
        name = "gibbs"
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("compute the g/(RT) at the given temperature")
    )
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    if syms is None:
        generate_thermo_routine(
            fstream,
            species_info,
            name,
            gibbs_nasa,
            species_coeffs,
            qss_flag,
            1,
        )
    else:
        generate_thermo_routine(
            fstream,
            species_info,
            name,
            gibbs_nasa,
            species_coeffs,
            qss_flag,
            1,
            syms,
        )


def helmholtz(fstream, species_info, species_coeffs):
    """Write Helmholtz."""
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("compute the a/(RT) at the given temperature")
    )
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    generate_thermo_routine(
        fstream,
        species_info,
        "helmholtz",
        helmholtz_nasa,
        species_coeffs,
        0,
        1,
    )


def species_internal_energy(fstream, species_info, species_coeffs):
    """Write species internal energy."""
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("compute the e/(RT) at the given temperature")
    )
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    generate_thermo_routine(
        fstream,
        species_info,
        "speciesInternalEnergy",
        internal_energy,
        species_coeffs,
        0,
        1,
    )


def species_enthalpy(
    fstream, species_info, species_coeffs, qss_flag, syms=None
):
    """Write species enthalpy."""
    if qss_flag:
        name = "speciesEnthalpy_qss"
    else:
        name = "speciesEnthalpy"
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("compute the h/(RT) at the given temperature (Eq 20)"),
    )
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )

    if syms is None:
        generate_thermo_routine(
            fstream,
            species_info,
            name,
            enthalpy_nasa,
            species_coeffs,
            qss_flag,
            1,
        )
    else:
        generate_thermo_routine(
            fstream,
            species_info,
            name,
            enthalpy_nasa,
            species_coeffs,
            qss_flag,
            1,
            syms,
        )


def species_entropy(fstream, species_info, species_coeffs):
    """Write species entropy."""
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("compute the S/R at the given temperature (Eq 21)")
    )
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    generate_thermo_routine(
        fstream,
        species_info,
        "speciesEntropy",
        entropy_nasa,
        species_coeffs,
        0,
    )


def dthermodtemp(fstream, mechanism, species_info):
    """Write gradient of thermodynamics quantity wrt temperature."""
    species_coeffs = analyze_thermodynamics(mechanism, species_info, 0)
    dcvpdtemp(fstream, species_info, species_coeffs)


def dcvpdtemp(fstream, species_info, species_coeffs):
    """Write gradient of cp/cv wrt temperature."""
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature"
        ),
    )
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    generate_thermo_routine(
        fstream, species_info, "dcvpRdT", dcpdtemp_nasa, species_coeffs, 0
    )


def dcpdtemp_nasa(fstream, parameters):
    """Write NASA polynomial for dcpdtemp."""
    cw.writer(fstream, "%+15.8e" % parameters[1])
    cw.writer(fstream, "%+15.8e * tc[1]" % (parameters[2] * 2.0))
    cw.writer(fstream, "%+15.8e * tc[2]" % (parameters[3] * 3.0))
    cw.writer(fstream, "%+15.8e * tc[3];" % (parameters[4] * 4.0))


def cv_nasa(fstream, parameters):
    """Write NASA polynomial for cv."""
    cw.writer(fstream, "%+15.8e" % (parameters[0] - 1.0))
    cw.writer(fstream, "%+15.8e * tc[1]" % parameters[1])
    cw.writer(fstream, "%+15.8e * tc[2]" % parameters[2])
    cw.writer(fstream, "%+15.8e * tc[3]" % parameters[3])
    cw.writer(fstream, "%+15.8e * tc[4];" % parameters[4])


def cp_nasa(fstream, parameters):
    """Write NASA polynomial for cp."""
    cw.writer(fstream, "%+15.8e" % parameters[0])
    cw.writer(fstream, "%+15.8e * tc[1]" % parameters[1])
    cw.writer(fstream, "%+15.8e * tc[2]" % parameters[2])
    cw.writer(fstream, "%+15.8e * tc[3]" % parameters[3])
    cw.writer(fstream, "%+15.8e * tc[4];" % parameters[4])


def gibbs_nasa(fstream, parameters, syms=None):
    """Write NASA polynomial for Gibbs."""
    record_symbolic_operations = True
    if syms is None:
        record_symbolic_operations = False

    if record_symbolic_operations:
        symb_smp = 0.0

    cw.writer(fstream, "%+20.15e * invT" % parameters[5])
    if record_symbolic_operations:
        symb_smp += parameters[5] * syms.invT_smp
    cw.writer(fstream, "%+20.15e" % (parameters[0] - parameters[6]))
    if record_symbolic_operations:
        symb_smp += parameters[0] - parameters[6]
    cw.writer(fstream, "%+20.15e * tc[0]" % (-parameters[0]))
    if record_symbolic_operations:
        symb_smp += (-parameters[0]) * syms.tc_smp[0]
    cw.writer(fstream, "%+20.15e * tc[1]" % ((-parameters[1] / 2)))
    if record_symbolic_operations:
        symb_smp += (-parameters[1] / 2) * syms.tc_smp[1]
    cw.writer(fstream, "%+20.15e * tc[2]" % ((-parameters[2] / 6)))
    if record_symbolic_operations:
        symb_smp += (-parameters[2] / 6) * syms.tc_smp[2]
    cw.writer(fstream, "%+20.15e * tc[3]" % ((-parameters[3] / 12)))
    if record_symbolic_operations:
        symb_smp += (-parameters[3] / 12) * syms.tc_smp[3]
    cw.writer(fstream, "%+20.15e * tc[4];" % ((-parameters[4] / 20)))
    if record_symbolic_operations:
        symb_smp += (-parameters[4] / 20) * syms.tc_smp[4]

    if record_symbolic_operations:
        return symb_smp


def helmholtz_nasa(fstream, parameters, syms=None):
    """Write NASA polynomial for Helmholtz."""
    record_symbolic_operations = True
    if syms is None:
        record_symbolic_operations = False

    if record_symbolic_operations:
        symb_smp = 0.0

    cw.writer(fstream, "%+15.8e * invT" % parameters[5])
    if record_symbolic_operations:
        symb_smp += parameters[5] * syms.invT_smp
    cw.writer(fstream, "%+15.8e" % (parameters[0] - parameters[6] - 1.0))
    if record_symbolic_operations:
        symb_smp += parameters[0] - parameters[6] - 1.0
    cw.writer(fstream, "%+15.8e * tc[0]" % (-parameters[0]))
    if record_symbolic_operations:
        symb_smp += (-parameters[0]) * syms.tc_smp[0]
    cw.writer(fstream, "%+15.8e * tc[1]" % ((-parameters[1] / 2)))
    if record_symbolic_operations:
        symb_smp += (-parameters[1] / 2) * syms.tc_smp[1]
    cw.writer(fstream, "%+15.8e * tc[2]" % ((-parameters[2] / 6)))
    if record_symbolic_operations:
        symb_smp += (-parameters[2] / 6) * syms.tc_smp[2]
    cw.writer(fstream, "%+15.8e * tc[3]" % ((-parameters[3] / 12)))
    if record_symbolic_operations:
        symb_smp += (-parameters[3] / 12) * syms.tc_smp[3]
    cw.writer(fstream, "%+15.8e * tc[4];" % ((-parameters[4] / 20)))
    if record_symbolic_operations:
        symb_smp += (-parameters[4] / 20) * syms.tc_smp[4]

    if record_symbolic_operations:
        return symb_smp


def internal_energy(fstream, parameters, syms=None):
    """Write NASA polynomial for internal energy."""
    record_symbolic_operations = True
    if syms is None:
        record_symbolic_operations = False

    if record_symbolic_operations:
        symb_smp = 0.0

    cw.writer(fstream, "%+15.8e" % (parameters[0] - 1.0))
    if record_symbolic_operations:
        symb_smp += parameters[0] - 1.0
    cw.writer(fstream, "%+15.8e * tc[1]" % ((parameters[1] / 2)))
    if record_symbolic_operations:
        symb_smp += (parameters[1] / 2) * syms.tc_smp[1]
    cw.writer(fstream, "%+15.8e * tc[2]" % ((parameters[2] / 3)))
    if record_symbolic_operations:
        symb_smp += (parameters[2] / 3) * syms.tc_smp[2]
    cw.writer(fstream, "%+15.8e * tc[3]" % ((parameters[3] / 4)))
    if record_symbolic_operations:
        symb_smp += (parameters[3] / 4) * syms.tc_smp[3]
    cw.writer(fstream, "%+15.8e * tc[4]" % ((parameters[4] / 5)))
    if record_symbolic_operations:
        symb_smp += (parameters[4] / 5) * syms.tc_smp[4]
    cw.writer(fstream, "%+15.8e * invT;" % (parameters[5]))
    if record_symbolic_operations:
        symb_smp += (parameters[5]) * syms.invT_smp

    if record_symbolic_operations:
        return symb_smp


def enthalpy_nasa(fstream, parameters, syms=None):
    """Write NASA polynomial for enthalpy."""
    record_symbolic_operations = True
    if syms is None:
        record_symbolic_operations = False

    if record_symbolic_operations:
        symb_smp = 0.0

    cw.writer(fstream, "%+15.8e" % parameters[0])
    if record_symbolic_operations:
        symb_smp += parameters[0]
    cw.writer(fstream, "%+15.8e * tc[1]" % ((parameters[1] / 2)))
    if record_symbolic_operations:
        symb_smp += (parameters[1] / 2) * syms.tc_smp[1]
    cw.writer(fstream, "%+15.8e * tc[2]" % ((parameters[2] / 3)))
    if record_symbolic_operations:
        symb_smp += (parameters[2] / 3) * syms.tc_smp[2]
    cw.writer(fstream, "%+15.8e * tc[3]" % ((parameters[3] / 4)))
    if record_symbolic_operations:
        symb_smp += (parameters[3] / 4) * syms.tc_smp[3]
    cw.writer(fstream, "%+15.8e * tc[4]" % ((parameters[4] / 5)))
    if record_symbolic_operations:
        symb_smp += (parameters[4] / 5) * syms.tc_smp[4]
    cw.writer(fstream, "%+15.8e * invT;" % (parameters[5]))
    if record_symbolic_operations:
        symb_smp += (parameters[5]) * syms.invT_smp

    if record_symbolic_operations:
        return symb_smp


def entropy_nasa(fstream, parameters, syms=None):
    """Write NASA polynomial for entropy."""
    record_symbolic_operations = True
    if syms is None:
        record_symbolic_operations = False

    if record_symbolic_operations:
        symb_smp = 0.0

    cw.writer(fstream, "%+15.8e * tc[0]" % parameters[0])
    if record_symbolic_operations:
        symb_smp += (parameters[0]) * syms.tc_smp[0]
    cw.writer(fstream, "%+15.8e * tc[1]" % (parameters[1]))
    if record_symbolic_operations:
        symb_smp += (parameters[1]) * syms.tc_smp[1]
    cw.writer(fstream, "%+15.8e * tc[2]" % ((parameters[2] / 2)))
    if record_symbolic_operations:
        symb_smp += (parameters[2] / 2) * syms.tc_smp[2]
    cw.writer(fstream, "%+15.8e * tc[3]" % ((parameters[3] / 3)))
    if record_symbolic_operations:
        symb_smp += (parameters[3] / 3) * syms.tc_smp[3]
    cw.writer(fstream, "%+15.8e * tc[4]" % ((parameters[4] / 4)))
    if record_symbolic_operations:
        symb_smp += (parameters[4] / 4) * syms.tc_smp[4]
    cw.writer(fstream, "%+15.8e ;" % (parameters[6]))
    if record_symbolic_operations:
        symb_smp += parameters[6]

    if record_symbolic_operations:
        return symb_smp
