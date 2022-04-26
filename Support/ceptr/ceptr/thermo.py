import sys
from collections import OrderedDict

import ceptr.writer as cw


def thermo(fstream, mechanism, species_info):
    species_coeffs = analyze_thermodynamics(mechanism, species_info, 0)
    if species_info.nQSSspecies > 0:
        qss_species_coeffs = analyze_thermodynamics(mechanism, species_info, 1)

    cv(fstream, species_info, species_coeffs)
    cp(fstream, species_info, species_coeffs)
    gibbs(fstream, species_info, species_coeffs, 0)
    if species_info.nQSSspecies > 0:
        gibbs(fstream, species_info, qss_species_coeffs, 1)
    helmholtz(fstream, species_info, species_coeffs)
    species_internal_energy(fstream, species_info, species_coeffs)
    species_enthalpy(fstream, species_info, species_coeffs, 0)
    if species_info.nQSSspecies > 0:
        species_enthalpy(fstream, species_info, qss_species_coeffs, 1)
    species_entropy(fstream, species_info, species_coeffs)


def analyze_thermodynamics(mechanism, species_info, qss_flag):
    low_temp = 0.0
    high_temp = 1000000.0

    midpoints = OrderedDict()

    if qss_flag:
        for symbol in species_info.qss_species_list:
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
        for symbol in species_info.nonqss_species_list:
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
):
    low_temp, high_temp, midpoints = species_coeffs

    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void %s(amrex::Real *"
        " species, const amrex::Real *  tc)" % name,
    )
    cw.writer(fstream, "{")
    # declarations
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("temperature"))
    cw.writer(fstream, "const amrex::Real T = tc[1];")
    if needs_inv_temp != 0:
        cw.writer(fstream, "const amrex::Real invT = 1 / T;")
    if needs_inv_temp == 2:
        cw.writer(fstream, "const amrex::Real invT2 = invT*invT;")

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
            expression_generator(fstream, high_range)

        cw.writer(fstream, "}")

    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def cv(fstream, species_info, species_coeffs):
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
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute Cp/R at the given temperature"))
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    generate_thermo_routine(
        fstream, species_info, "cp_R", cp_nasa, species_coeffs, 0
    )


def gibbs(fstream, species_info, species_coeffs, qss_flag):
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
    generate_thermo_routine(
        fstream, species_info, name, gibbs_nasa, species_coeffs, qss_flag, 1
    )


def helmholtz(fstream, species_info, species_coeffs):
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


def species_enthalpy(fstream, species_info, species_coeffs, qss_flag):
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
    generate_thermo_routine(
        fstream, species_info, name, enthalpy_nasa, species_coeffs, qss_flag, 1
    )


def species_entropy(fstream, species_info, species_coeffs):
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
    species_coeffs = analyze_thermodynamics(mechanism, species_info, 0)
    dcvpdtemp(fstream, species_info, species_coeffs)


def dcvpdtemp(fstream, species_info, species_coeffs):
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
    cw.writer(fstream, "%+15.8e" % parameters[1])
    cw.writer(fstream, "%+15.8e * tc[1]" % (parameters[2] * 2.0))
    cw.writer(fstream, "%+15.8e * tc[2]" % (parameters[3] * 3.0))
    cw.writer(fstream, "%+15.8e * tc[3];" % (parameters[4] * 4.0))


def cv_nasa(fstream, parameters):
    cw.writer(fstream, "%+15.8e" % (parameters[0] - 1.0))
    cw.writer(fstream, "%+15.8e * tc[1]" % parameters[1])
    cw.writer(fstream, "%+15.8e * tc[2]" % parameters[2])
    cw.writer(fstream, "%+15.8e * tc[3]" % parameters[3])
    cw.writer(fstream, "%+15.8e * tc[4];" % parameters[4])


def cp_nasa(fstream, parameters):
    cw.writer(fstream, "%+15.8e" % parameters[0])
    cw.writer(fstream, "%+15.8e * tc[1]" % parameters[1])
    cw.writer(fstream, "%+15.8e * tc[2]" % parameters[2])
    cw.writer(fstream, "%+15.8e * tc[3]" % parameters[3])
    cw.writer(fstream, "%+15.8e * tc[4];" % parameters[4])


def gibbs_nasa(fstream, parameters):
    cw.writer(fstream, "%+20.15e * invT" % parameters[5])
    cw.writer(fstream, "%+20.15e" % (parameters[0] - parameters[6]))
    cw.writer(fstream, "%+20.15e * tc[0]" % (-parameters[0]))
    cw.writer(fstream, "%+20.15e * tc[1]" % ((-parameters[1] / 2)))
    cw.writer(fstream, "%+20.15e * tc[2]" % ((-parameters[2] / 6)))
    cw.writer(fstream, "%+20.15e * tc[3]" % ((-parameters[3] / 12)))
    cw.writer(fstream, "%+20.15e * tc[4];" % ((-parameters[4] / 20)))


def helmholtz_nasa(fstream, parameters):
    cw.writer(fstream, "%+15.8e * invT" % parameters[5])
    cw.writer(fstream, "%+15.8e" % (parameters[0] - parameters[6] - 1.0))
    cw.writer(fstream, "%+15.8e * tc[0]" % (-parameters[0]))
    cw.writer(fstream, "%+15.8e * tc[1]" % ((-parameters[1] / 2)))
    cw.writer(fstream, "%+15.8e * tc[2]" % ((-parameters[2] / 6)))
    cw.writer(fstream, "%+15.8e * tc[3]" % ((-parameters[3] / 12)))
    cw.writer(fstream, "%+15.8e * tc[4];" % ((-parameters[4] / 20)))


def internal_energy(fstream, parameters):
    cw.writer(fstream, "%+15.8e" % (parameters[0] - 1.0))
    cw.writer(fstream, "%+15.8e * tc[1]" % ((parameters[1] / 2)))
    cw.writer(fstream, "%+15.8e * tc[2]" % ((parameters[2] / 3)))
    cw.writer(fstream, "%+15.8e * tc[3]" % ((parameters[3] / 4)))
    cw.writer(fstream, "%+15.8e * tc[4]" % ((parameters[4] / 5)))
    cw.writer(fstream, "%+15.8e * invT;" % (parameters[5]))


def enthalpy_nasa(fstream, parameters):
    cw.writer(fstream, "%+15.8e" % parameters[0])
    cw.writer(fstream, "%+15.8e * tc[1]" % ((parameters[1] / 2)))
    cw.writer(fstream, "%+15.8e * tc[2]" % ((parameters[2] / 3)))
    cw.writer(fstream, "%+15.8e * tc[3]" % ((parameters[3] / 4)))
    cw.writer(fstream, "%+15.8e * tc[4]" % ((parameters[4] / 5)))
    cw.writer(fstream, "%+15.8e * invT;" % (parameters[5]))


def entropy_nasa(fstream, parameters):
    cw.writer(fstream, "%+15.8e * tc[0]" % parameters[0])
    cw.writer(fstream, "%+15.8e * tc[1]" % (parameters[1]))
    cw.writer(fstream, "%+15.8e * tc[2]" % ((parameters[2] / 2)))
    cw.writer(fstream, "%+15.8e * tc[3]" % ((parameters[3] / 3)))
    cw.writer(fstream, "%+15.8e * tc[4]" % ((parameters[4] / 4)))
    cw.writer(fstream, "%+15.8e ;" % (parameters[6]))
