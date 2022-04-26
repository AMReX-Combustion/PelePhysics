import sys
from collections import OrderedDict

import ceptr.writer as cw


def thermo(fstream, mechanism, species_info):
    species_coeffs = analyzeThermodynamics(mechanism, species_info, 0)
    if species_info.nQSSspecies > 0:
        QSSspecies_coeffs = analyzeThermodynamics(mechanism, species_info, 1)

    cv(fstream, species_info, species_coeffs)
    cp(fstream, species_info, species_coeffs)
    gibbs(fstream, species_info, species_coeffs, 0)
    if species_info.nQSSspecies > 0:
        gibbs(fstream, species_info, QSSspecies_coeffs, 1)
    helmholtz(fstream, species_info, species_coeffs)
    species_internal_energy(fstream, species_info, species_coeffs)
    speciesEnthalpy(fstream, species_info, species_coeffs, 0)
    if species_info.nQSSspecies > 0:
        speciesEnthalpy(fstream, species_info, QSSspecies_coeffs, 1)
    species_entropy(fstream, species_info, species_coeffs)


def analyzeThermodynamics(mechanism, species_info, qss_flag):
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


def generateThermoRoutine(
    fstream,
    species_info,
    name,
    expressionGenerator,
    species_coeffs,
    QSS_flag,
    needsInvT=0,
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
    if needsInvT != 0:
        cw.writer(fstream, "const amrex::Real invT = 1 / T;")
    if needsInvT == 2:
        cw.writer(fstream, "const amrex::Real invT2 = invT*invT;")

    for mid_temp, speciesList in list(midpoints.items()):
        cw.writer(fstream, "")
        cw.writer(
            fstream,
            cw.comment("species with midpoint at T=%g kelvin" % mid_temp),
        )
        cw.writer(fstream, "if (T < %g) {" % mid_temp)

        for species, low_range, _ in speciesList:
            if QSS_flag:
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
            expressionGenerator(fstream, low_range)

        cw.writer(fstream, "} else {")

        for species, _, high_range in speciesList:
            if QSS_flag:
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
            expressionGenerator(fstream, high_range)

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
    generateThermoRoutine(
        fstream, species_info, "cv_R", cvNASA, species_coeffs, 0
    )


def cp(fstream, species_info, species_coeffs):
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute Cp/R at the given temperature"))
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    generateThermoRoutine(
        fstream, species_info, "cp_R", cpNASA, species_coeffs, 0
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
    generateThermoRoutine(
        fstream, species_info, name, gibbsNASA, species_coeffs, qss_flag, 1
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
    generateThermoRoutine(
        fstream, species_info, "helmholtz", helmholtzNASA, species_coeffs, 0, 1
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
    generateThermoRoutine(
        fstream,
        species_info,
        "speciesInternalEnergy",
        internalEnergy,
        species_coeffs,
        0,
        1,
    )


def speciesEnthalpy(fstream, species_info, species_coeffs, qss_flag):
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
    generateThermoRoutine(
        fstream, species_info, name, enthalpyNASA, species_coeffs, qss_flag, 1
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
    generateThermoRoutine(
        fstream, species_info, "speciesEntropy", entropyNASA, species_coeffs, 0
    )


def dthermodT(fstream, mechanism, species_info):
    species_coeffs = analyzeThermodynamics(mechanism, species_info, 0)
    dcvpdT(fstream, species_info, species_coeffs)


def dcvpdT(fstream, species_info, species_coeffs):
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
    generateThermoRoutine(
        fstream, species_info, "dcvpRdT", dcpdTNASA, species_coeffs, 0
    )


def dcpdTNASA(fstream, parameters):
    cw.writer(fstream, "%+15.8e" % parameters[1])
    cw.writer(fstream, "%+15.8e * tc[1]" % (parameters[2] * 2.0))
    cw.writer(fstream, "%+15.8e * tc[2]" % (parameters[3] * 3.0))
    cw.writer(fstream, "%+15.8e * tc[3];" % (parameters[4] * 4.0))


def cvNASA(fstream, parameters):
    cw.writer(fstream, "%+15.8e" % (parameters[0] - 1.0))
    cw.writer(fstream, "%+15.8e * tc[1]" % parameters[1])
    cw.writer(fstream, "%+15.8e * tc[2]" % parameters[2])
    cw.writer(fstream, "%+15.8e * tc[3]" % parameters[3])
    cw.writer(fstream, "%+15.8e * tc[4];" % parameters[4])


def cpNASA(fstream, parameters):
    cw.writer(fstream, "%+15.8e" % parameters[0])
    cw.writer(fstream, "%+15.8e * tc[1]" % parameters[1])
    cw.writer(fstream, "%+15.8e * tc[2]" % parameters[2])
    cw.writer(fstream, "%+15.8e * tc[3]" % parameters[3])
    cw.writer(fstream, "%+15.8e * tc[4];" % parameters[4])


def gibbsNASA(fstream, parameters):
    cw.writer(fstream, "%+20.15e * invT" % parameters[5])
    cw.writer(fstream, "%+20.15e" % (parameters[0] - parameters[6]))
    cw.writer(fstream, "%+20.15e * tc[0]" % (-parameters[0]))
    cw.writer(fstream, "%+20.15e * tc[1]" % ((-parameters[1] / 2)))
    cw.writer(fstream, "%+20.15e * tc[2]" % ((-parameters[2] / 6)))
    cw.writer(fstream, "%+20.15e * tc[3]" % ((-parameters[3] / 12)))
    cw.writer(fstream, "%+20.15e * tc[4];" % ((-parameters[4] / 20)))


def helmholtzNASA(fstream, parameters):
    cw.writer(fstream, "%+15.8e * invT" % parameters[5])
    cw.writer(fstream, "%+15.8e" % (parameters[0] - parameters[6] - 1.0))
    cw.writer(fstream, "%+15.8e * tc[0]" % (-parameters[0]))
    cw.writer(fstream, "%+15.8e * tc[1]" % ((-parameters[1] / 2)))
    cw.writer(fstream, "%+15.8e * tc[2]" % ((-parameters[2] / 6)))
    cw.writer(fstream, "%+15.8e * tc[3]" % ((-parameters[3] / 12)))
    cw.writer(fstream, "%+15.8e * tc[4];" % ((-parameters[4] / 20)))


def internalEnergy(fstream, parameters):
    cw.writer(fstream, "%+15.8e" % (parameters[0] - 1.0))
    cw.writer(fstream, "%+15.8e * tc[1]" % ((parameters[1] / 2)))
    cw.writer(fstream, "%+15.8e * tc[2]" % ((parameters[2] / 3)))
    cw.writer(fstream, "%+15.8e * tc[3]" % ((parameters[3] / 4)))
    cw.writer(fstream, "%+15.8e * tc[4]" % ((parameters[4] / 5)))
    cw.writer(fstream, "%+15.8e * invT;" % (parameters[5]))


def enthalpyNASA(fstream, parameters):
    cw.writer(fstream, "%+15.8e" % parameters[0])
    cw.writer(fstream, "%+15.8e * tc[1]" % ((parameters[1] / 2)))
    cw.writer(fstream, "%+15.8e * tc[2]" % ((parameters[2] / 3)))
    cw.writer(fstream, "%+15.8e * tc[3]" % ((parameters[3] / 4)))
    cw.writer(fstream, "%+15.8e * tc[4]" % ((parameters[4] / 5)))
    cw.writer(fstream, "%+15.8e * invT;" % (parameters[5]))


def entropyNASA(fstream, parameters):
    cw.writer(fstream, "%+15.8e * tc[0]" % parameters[0])
    cw.writer(fstream, "%+15.8e * tc[1]" % (parameters[1]))
    cw.writer(fstream, "%+15.8e * tc[2]" % ((parameters[2] / 2)))
    cw.writer(fstream, "%+15.8e * tc[3]" % ((parameters[3] / 3)))
    cw.writer(fstream, "%+15.8e * tc[4]" % ((parameters[4] / 4)))
    cw.writer(fstream, "%+15.8e ;" % (parameters[6]))
