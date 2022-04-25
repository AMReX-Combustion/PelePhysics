import sys
from collections import OrderedDict

import ceptr.writer as cw


def thermo(fstream, mechanism, species_info):
    speciesCoeffs = analyzeThermodynamics(mechanism, species_info, 0)
    if species_info.nQSSspecies > 0:
        QSSspeciesCoeffs = analyzeThermodynamics(mechanism, species_info, 1)

    cv(fstream, species_info, speciesCoeffs)
    cp(fstream, species_info, speciesCoeffs)
    gibbs(fstream, species_info, speciesCoeffs, 0)
    if species_info.nQSSspecies > 0:
        gibbs(fstream, species_info, QSSspeciesCoeffs, 1)
    helmholtz(fstream, species_info, speciesCoeffs)
    speciesInternalEnergy(fstream, species_info, speciesCoeffs)
    speciesEnthalpy(fstream, species_info, speciesCoeffs, 0)
    if species_info.nQSSspecies > 0:
        speciesEnthalpy(fstream, species_info, QSSspeciesCoeffs, 1)
    speciesEntropy(fstream, species_info, speciesCoeffs)


def analyzeThermodynamics(mechanism, species_info, QSS_Flag):
    lowT = 0.0
    highT = 1000000.0

    midpoints = OrderedDict()

    lowT_qss = 0.0
    highT_qss = 1000000.0

    midpoints_qss = {}

    # FIXME
    if QSS_Flag:
        for symbol in self.qss_species_list:
            species = self.mechanism.species(symbol)
            models = species.thermo
            if len(models) > 2:
                print("species: ", species)
                import pyre

                pyre.debug.Firewall.hit(
                    "unsupported configuration in species.thermo"
                )
                return

            m1 = models[0]
            m2 = models[1]

            if m1.lowT < m2.lowT:
                lowRange = m1
                highRange = m2
            else:
                lowRange = m2
                highRange = m1

            low = lowRange.lowT
            mid = lowRange.highT
            high = highRange.highT

            if low > lowT:
                lowT = low
            if high < highT:
                highT = high

            midpoints.setdefault(mid, []).append(
                (species, lowRange, highRange)
            )

    else:
        for symbol in species_info.nonqss_species_list:
            species = mechanism.species(symbol)
            model = species.thermo

            if not model.n_coeffs == 15:
                print("Unsupported thermo model.")
                sys.exit(1)

            loT = model.min_temp
            hiT = model.max_temp
            if lowT < loT:
                lowT = loT
            if hiT < highT:
                highT = hiT
            mid = model.coeffs[0]
            highRange = model.coeffs[1:8]
            lowRange = model.coeffs[8:15]

            midpoints.setdefault(mid, []).append(
                (species, lowRange, highRange)
            )

    species_info.lowT = lowT
    species_info.highT = highT
    return lowT, highT, midpoints


def generateThermoRoutine(
    fstream,
    species_info,
    name,
    expressionGenerator,
    speciesCoeffs,
    QSS_flag,
    needsInvT=0,
):
    lowT, highT, midpoints = speciesCoeffs

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

    for midT, speciesList in list(midpoints.items()):
        cw.writer(fstream, "")
        cw.writer(
            fstream, cw.comment("species with midpoint at T=%g kelvin" % midT)
        )
        cw.writer(fstream, "if (T < %g) {" % midT)

        for species, lowRange, highRange in speciesList:
            if QSS_flag:
                cw.writer(
                    fstream,
                    cw.comment(
                        "species %d: %s"
                        % (
                            species_info.ordered_idx_map[species.name]
                            - species_info.nSpecies,
                            species.name,
                        )
                    ),
                )
                cw.writer(
                    fstream,
                    "species[%d] ="
                    % (
                        species_info.ordered_idx_map[species.name]
                        - species_info.nSpecies
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
            expressionGenerator(fstream, lowRange)

        cw.writer(fstream, "} else {")

        for species, lowRange, highRange in speciesList:
            if QSS_flag:
                cw.writer(
                    fstream,
                    cw.comment(
                        "species %d: %s"
                        % (
                            species_info.ordered_idx_map[species.name]
                            - species_info.nSpecies,
                            species.name,
                        )
                    ),
                )
                cw.writer(
                    fstream,
                    "species[%d] ="
                    % (
                        species_info.ordered_idx_map[species.name]
                        - species_info.nSpecies
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
            expressionGenerator(fstream, highRange)

        cw.writer(fstream, "}")

    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def cv(fstream, species_info, speciesCoeffs):
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute Cv/R at the given temperature"))
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    generateThermoRoutine(
        fstream, species_info, "cv_R", cvNASA, speciesCoeffs, 0
    )


def cp(fstream, species_info, speciesCoeffs):
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute Cp/R at the given temperature"))
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    generateThermoRoutine(
        fstream, species_info, "cp_R", cpNASA, speciesCoeffs, 0
    )


def gibbs(fstream, species_info, speciesCoeffs, QSS_Flag):
    if QSS_Flag:
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
        fstream, species_info, name, gibbsNASA, speciesCoeffs, QSS_Flag, 1
    )


def helmholtz(fstream, species_info, speciesCoeffs):
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("compute the a/(RT) at the given temperature")
    )
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    generateThermoRoutine(
        fstream, species_info, "helmholtz", helmholtzNASA, speciesCoeffs, 0, 1
    )


def speciesInternalEnergy(fstream, species_info, speciesCoeffs):
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
        speciesCoeffs,
        0,
        1,
    )


def speciesEnthalpy(fstream, species_info, speciesCoeffs, QSS_Flag):
    if QSS_Flag:
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
        fstream, species_info, name, enthalpyNASA, speciesCoeffs, QSS_Flag, 1
    )


def speciesEntropy(fstream, species_info, speciesCoeffs):
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("compute the S/R at the given temperature (Eq 21)")
    )
    cw.writer(
        fstream,
        cw.comment("tc contains precomputed powers of T, tc[0] = log(T)"),
    )
    generateThermoRoutine(
        fstream, species_info, "speciesEntropy", entropyNASA, speciesCoeffs, 0
    )


def dthermodT(fstream, mechanism, species_info):
    speciesCoeffs = analyzeThermodynamics(mechanism, species_info, 0)
    dcvpdT(fstream, species_info, speciesCoeffs)


def dcvpdT(fstream, species_info, speciesCoeffs):
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
        fstream, species_info, "dcvpRdT", dcpdTNASA, speciesCoeffs, 0
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
