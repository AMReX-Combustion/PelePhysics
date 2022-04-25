import sys
from collections import OrderedDict

import numpy as np

import ceptr.constants as cc
import ceptr.writer as cw


def transport(fstream, mechanism, species_info):
    """Write the transport functions."""
    cw.writer(fstream, cw.comment("Transport function declarations "))
    nSpecies = species_info.nSpecies
    speciesTransport = analyzeTransport(mechanism, species_info)
    NLITE = 0
    idxLightSpecs = []
    for sp in range(nSpecies):
        spec = species_info.nonqss_species[sp]
        if spec.weight < 5.0:
            NLITE += 1
            idxLightSpecs.append(spec.idx)
    miscTransInfo(fstream, KK=nSpecies, NLITE=NLITE, do_declarations=False)
    wt(fstream, species_info, False)
    eps(fstream, mechanism, species_info, speciesTransport, False)
    sig(fstream, mechanism, species_info, speciesTransport, False)
    dip(fstream, mechanism, species_info, speciesTransport, False)
    pol(fstream, mechanism, species_info, speciesTransport, False)
    zrot(fstream, mechanism, species_info, speciesTransport, False)
    nlin(fstream, mechanism, species_info, speciesTransport, False)

    viscosity(
        fstream, mechanism, species_info, speciesTransport, False, NTFit=50
    )
    diffcoefs(fstream, species_info, speciesTransport, False, NTFit=50)
    lightSpecs(fstream, idxLightSpecs, False)
    thermaldiffratios(
        fstream, species_info, speciesTransport, idxLightSpecs, False, NTFit=50
    )


def analyzeTransport(mechanism, species_info):
    """Extract transport model coefficients."""
    transdata = OrderedDict()

    for spec in species_info.nonqss_species:
        species = mechanism.species(spec.name)

        m1 = species.transport
        if m1.geometry == "atom":
            lin = 0
        elif m1.geometry == "linear":
            lin = 1
        elif m1.geometry == "nonlinear":
            lin = 2
        else:
            print("Unrecognized species geometry in transport")
            sys.exit(1)
        eps = (m1.well_depth * cc.ureg.joule).to(cc.ureg.erg).m / cc.kb
        sig = (m1.diameter * cc.ureg.meter).to(cc.ureg.angstrom).m
        dip = (m1.dipole * cc.ureg.coulomb * cc.ureg.m).to(cc.ureg.debye).m
        pol = (
            (m1.polarizability * cc.ureg.meter**3)
            .to(cc.ureg.angstrom**3)
            .m
        )
        zrot = m1.rotational_relaxation

        transdata[spec] = [lin, eps, sig, dip, pol, zrot]

    return transdata


def miscTransInfo(fstream, KK, NLITE, do_declarations, NO=4):
    """Write transport information."""
    cw.writer(fstream)
    LENIMC = 4 * KK + NLITE
    generateTransRoutineInteger(
        fstream,
        [
            "egtransetLENIMC",
            "EGTRANSETLENIMC",
            "egtransetlenimc",
            "egtransetlenimc_",
            "LENIMC",
        ],
        LENIMC,
        do_declarations,
    )

    cw.writer(fstream)
    cw.writer(fstream)
    LENRMC = (19 + 2 * NO + NO * NLITE) * KK + (15 + NO) * KK**2
    generateTransRoutineInteger(
        fstream,
        [
            "egtransetLENRMC",
            "EGTRANSETLENRMC",
            "egtransetlenrmc",
            "egtransetlenrmc_",
            "LENRMC",
        ],
        LENRMC,
        do_declarations,
    )

    cw.writer(fstream)
    cw.writer(fstream)
    generateTransRoutineInteger(
        fstream,
        [
            "egtransetNO",
            "EGTRANSETNO",
            "egtransetno",
            "egtransetno_",
            "NO",
        ],
        NO,
        do_declarations,
    )

    cw.writer(fstream)
    cw.writer(fstream)
    generateTransRoutineInteger(
        fstream,
        [
            "egtransetKK",
            "EGTRANSETKK",
            "egtransetkk",
            "egtransetkk_",
            "KK",
        ],
        KK,
        do_declarations,
    )

    cw.writer(fstream)
    cw.writer(fstream)
    generateTransRoutineInteger(
        fstream,
        [
            "egtransetNLITE",
            "EGTRANSETNLITE",
            "egtransetnlite",
            "egtransetnlite_",
            "NLITE",
        ],
        NLITE,
        do_declarations,
    )

    cw.writer(fstream)
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Patm in ergs/cm3"))

    if do_declarations:
        cw.writer(fstream, "#if defined(BL_FORT_USE_UPPERCASE)")
        cw.writer(fstream, "#define egtransetPATM EGTRANSETPATM")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_LOWERCASE)")
        cw.writer(fstream, "#define egtransetPATM egtransetpatm")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_UNDERSCORE)")
        cw.writer(fstream, "#define egtransetPATM egtransetpatm_")
        cw.writer(fstream, "#endif")

    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
    cw.writer(fstream, "void egtransetPATM(amrex::Real* PATM) {")
    cw.writer(fstream, "*PATM =   0.1013250000000000E+07;}")


def generateTransRoutineInteger(fstream, nametab, expression, do_declarations):
    """Write generic integer transport routine."""
    if do_declarations:
        cw.writer(fstream, "#if defined(BL_FORT_USE_UPPERCASE)")
        cw.writer(fstream, "#define %s %s" % (nametab[0], nametab[1]))
        cw.writer(fstream, "#elif defined(BL_FORT_USE_LOWERCASE)")
        cw.writer(fstream, "#define %s %s" % (nametab[0], nametab[2]))
        cw.writer(fstream, "#elif defined(BL_FORT_USE_UNDERSCORE)")
        cw.writer(fstream, "#define %s %s" % (nametab[0], nametab[3]))
        cw.writer(fstream, "#endif")

    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
    cw.writer(fstream, "void %s(int* %s ) {" % (nametab[0], nametab[4]))

    cw.writer(fstream, "*%s = %d;}" % (nametab[4], expression))


def generateTransRoutineSimple(
    fstream,
    mechanism,
    species_info,
    nametab,
    idx,
    speciesTransport,
    do_declarations,
):
    """Write generic transport routine."""
    if do_declarations:
        cw.writer(fstream, "#if defined(BL_FORT_USE_UPPERCASE)")
        cw.writer(fstream, "#define %s %s" % (nametab[0], nametab[1]))
        cw.writer(fstream, "#elif defined(BL_FORT_USE_LOWERCASE)")
        cw.writer(fstream, "#define %s %s" % (nametab[0], nametab[2]))
        cw.writer(fstream, "#elif defined(BL_FORT_USE_UNDERSCORE)")
        cw.writer(fstream, "#define %s %s" % (nametab[0], nametab[3]))
        cw.writer(fstream, "#endif")

    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
    cw.writer(
        fstream, "void %s(amrex::Real* %s ) {" % (nametab[0], nametab[4])
    )

    for spec in species_info.nonqss_species:
        cw.writer(
            fstream,
            "%s[%d] = %.8E;"
            % (nametab[4], spec.idx, float(speciesTransport[spec][idx])),
        )
    cw.writer(fstream, "}")


def wt(fstream, species_info, do_declarations):
    """Write molecular weights function."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("the molecular weights in g/mol"))

    if do_declarations:
        cw.writer(fstream, "#if defined(BL_FORT_USE_UPPERCASE)")
        cw.writer(fstream, "#define egtransetWT EGTRANSETWT")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_LOWERCASE)")
        cw.writer(fstream, "#define egtransetWT egtransetwt")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_UNDERSCORE)")
        cw.writer(fstream, "#define egtransetWT egtransetwt_")
        cw.writer(fstream, "#endif")

    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
    cw.writer(fstream, "void %s(amrex::Real* %s ) {" % ("egtransetWT", "WT"))

    nSpecies = species_info.nSpecies
    for sp in range(nSpecies):
        species = species_info.nonqss_species[sp]
        cw.writer(
            fstream,
            "%s[%d] = %.8E;" % ("WT", species.idx, float(species.weight)),
        )

    cw.writer(fstream, "}")


def eps(fstream, mechanism, species_info, speciesTransport, do_declarations):
    """Write the lennard-jones potential well depth function."""
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("the lennard-jones potential well depth eps/kb in K"),
    )
    generateTransRoutineSimple(
        fstream,
        mechanism,
        species_info,
        [
            "egtransetEPS",
            "EGTRANSETEPS",
            "egtranseteps",
            "egtranseteps_",
            "EPS",
        ],
        1,
        speciesTransport,
        do_declarations,
    )


def sig(fstream, mechanism, species_info, speciesTransport, do_declarations):
    """Write the the lennard-jones collision diameter function."""
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("the lennard-jones collision diameter in Angstroms"),
    )
    generateTransRoutineSimple(
        fstream,
        mechanism,
        species_info,
        [
            "egtransetSIG",
            "EGTRANSETSIG",
            "egtransetsig",
            "egtransetsig_",
            "SIG",
        ],
        2,
        speciesTransport,
        do_declarations,
    )


def dip(fstream, mechanism, species_info, speciesTransport, do_declarations):
    """Write the dipole moment function."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("the dipole moment in Debye"))
    generateTransRoutineSimple(
        fstream,
        mechanism,
        species_info,
        [
            "egtransetDIP",
            "EGTRANSETDIP",
            "egtransetdip",
            "egtransetdip_",
            "DIP",
        ],
        3,
        speciesTransport,
        do_declarations,
    )


def pol(fstream, mechanism, species_info, speciesTransport, do_declarations):
    """Write the polarizability function."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("the polarizability in cubic Angstroms"))
    generateTransRoutineSimple(
        fstream,
        mechanism,
        species_info,
        [
            "egtransetPOL",
            "EGTRANSETPOL",
            "egtransetpol",
            "egtransetpol_",
            "POL",
        ],
        4,
        speciesTransport,
        do_declarations,
    )


def zrot(fstream, mechanism, species_info, speciesTransport, do_declarations):
    """Write the rotational relaxation collision number."""
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("the rotational relaxation collision number at 298 K"),
    )
    generateTransRoutineSimple(
        fstream,
        mechanism,
        species_info,
        [
            "egtransetZROT",
            "EGTRANSETZROT",
            "egtransetzrot",
            "egtransetzrot_",
            "ZROT",
        ],
        5,
        speciesTransport,
        do_declarations,
    )


def nlin(fstream, mechanism, species_info, speciesTransport, do_declarations):
    """Write the (monoatomic, linear, nonlinear) information."""
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("0: monoatomic, 1: linear, 2: nonlinear"))

    if do_declarations:
        cw.writer(fstream, "#if defined(BL_FORT_USE_UPPERCASE)")
        cw.writer(fstream, "#define egtransetNLIN EGTRANSETNLIN")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_LOWERCASE)")
        cw.writer(fstream, "#define egtransetNLIN egtransetnlin")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_UNDERSCORE)")
        cw.writer(fstream, "#define egtransetNLIN egtransetnlin_")
        cw.writer(fstream, "#endif")

    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
    cw.writer(fstream, "void egtransetNLIN(int* NLIN) {")

    for species in species_info.nonqss_species:
        cw.writer(
            fstream,
            "%s[%d] = %d;"
            % ("NLIN", species.idx, int(speciesTransport[species][0])),
        )

    cw.writer(fstream, "}")


def viscosity(
    fstream, mechanism, species_info, speciesTransport, do_declarations, NTFit
):
    """Write the viscosity function."""
    nSpecies = species_info.nSpecies
    # compute single constants in g/cm/s
    Na = 6.02214199e23
    RU = 8.31447e7
    # conversion coefs
    AtoCM = 1.0e-8
    DEBYEtoCGS = 1.0e-18
    # temperature increment
    dt = (species_info.highT - species_info.lowT) / (NTFit - 1)
    # factor dependent upon the molecule
    m_crot = np.zeros(nSpecies)
    m_cvib = np.zeros(nSpecies)
    isatm = np.zeros(nSpecies)
    for spec in speciesTransport:
        if int(speciesTransport[spec][0]) == 0:
            m_crot[spec.idx] = 0.0
            m_cvib[spec.idx] = 0.0
            isatm[spec.idx] = 0.0
        elif int(speciesTransport[spec][0]) == 1:
            m_crot[spec.idx] = 1.0
            m_cvib[spec.idx] = 5.0 / 2.0
            isatm[spec.idx] = 1.0
        else:
            m_crot[spec.idx] = 1.5
            m_cvib[spec.idx] = 3.0
            isatm[spec.idx] = 1.0
    # viscosities coefs (4 per spec)
    cofeta = OrderedDict()
    # conductivities coefs (4 per spec)
    coflam = OrderedDict()
    for spec in speciesTransport:
        spvisc = []
        spcond = []
        tlog = []
        for n in range(NTFit):
            t = species_info.lowT + dt * n
            # variables
            # eq. (2)
            tr = t / float(speciesTransport[spec][1])
            conversion = (
                DEBYEtoCGS * DEBYEtoCGS / AtoCM / AtoCM / AtoCM / cc.kb
            )
            dst = (
                0.5
                * conversion
                * float(speciesTransport[spec][3]) ** 2
                / (
                    float(speciesTransport[spec][1])
                    * float(speciesTransport[spec][2]) ** 3
                )
            )
            # viscosity of spec at t
            # eq. (1)
            conversion = AtoCM * AtoCM
            visc = (
                (5.0 / 16.0)
                * np.sqrt(np.pi * spec.weight * cc.kb * t / Na)
                / (
                    om22_CHEMKIN(tr, dst)
                    * np.pi
                    * float(speciesTransport[spec][2])
                    * float(speciesTransport[spec][2])
                    * conversion
                )
            )
            # conductivity of spec at t
            # eq. (30)
            conversion = AtoCM * AtoCM
            m_red = spec.weight / (2.0 * Na)
            diffcoef = (
                (3.0 / 16.0)
                * np.sqrt(2.0 * np.pi * cc.kb**3 * t**3 / m_red)
                / (
                    10.0
                    * np.pi
                    * om11_CHEMKIN(tr, dst)
                    * float(speciesTransport[spec][2])
                    * float(speciesTransport[spec][2])
                    * conversion
                )
            )
            # eq. (19)
            cv_vib_R = (
                getCVdRspecies(mechanism, t, spec) - m_cvib[spec.idx]
            ) * isatm[spec.idx]
            rho_atm = 10.0 * spec.weight / (RU * t)
            f_vib = rho_atm * diffcoef / visc
            # eq. (20)
            A = 2.5 - f_vib
            # eqs. (21) + (32-33)
            cv_rot_R = m_crot[spec.idx]
            # note: the T corr is not applied in CANTERA
            B = float(speciesTransport[spec][5]) * Fcorr(
                298.0, float(speciesTransport[spec][1])
            ) / Fcorr(t, float(speciesTransport[spec][1])) + (2.0 / np.pi) * (
                (5.0 / 3.0) * cv_rot_R + f_vib
            )
            # eq. (18)
            f_rot = f_vib * (1.0 + 2.0 / np.pi * A / B)
            # eq. (17)
            cv_trans_R = 3.0 / 2.0
            f_trans = (
                5.0 / 2.0 * (1.0 - 2.0 / np.pi * A / B * cv_rot_R / cv_trans_R)
            )
            if int(speciesTransport[spec][0]) == 0:
                cond = ((visc * RU / spec.weight)) * (5.0 / 2.0) * cv_trans_R
            else:
                cond = ((visc * RU / spec.weight)) * (
                    f_trans * cv_trans_R + f_rot * cv_rot_R + f_vib * cv_vib_R
                )

            # log transformation for polyfit
            tlog.append(np.log(t))
            spvisc.append(np.log(visc))
            spcond.append(np.log(cond))

        cofeta[spec.idx] = np.polyfit(tlog, spvisc, 3)
        coflam[spec.idx] = np.polyfit(tlog, spcond, 3)

    # header for visco
    cw.writer(fstream)
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Poly fits for the viscosities, dim NO*KK"))
    if do_declarations:
        cw.writer(fstream, "#if defined(BL_FORT_USE_UPPERCASE)")
        cw.writer(fstream, "#define egtransetCOFETA EGTRANSETCOFETA")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_LOWERCASE)")
        cw.writer(fstream, "#define egtransetCOFETA egtransetcofeta")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_UNDERSCORE)")
        cw.writer(fstream, "#define egtransetCOFETA egtransetcofeta_")
        cw.writer(fstream, "#endif")

    # visco coefs
    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
    cw.writer(fstream, "void egtransetCOFETA(amrex::Real* COFETA) {")

    for spec in species_info.nonqss_species:
        for i in range(4):
            cw.writer(
                fstream,
                "%s[%d] = %.8E;"
                % ("COFETA", spec.idx * 4 + i, cofeta[spec.idx][3 - i]),
            )

    cw.writer(fstream, "}")

    # header for cond
    cw.writer(fstream)
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Poly fits for the conductivities, dim NO*KK")
    )
    if do_declarations:
        cw.writer(fstream, "#if defined(BL_FORT_USE_UPPERCASE)")
        cw.writer(fstream, "#define egtransetCOFLAM EGTRANSETCOFLAM")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_LOWERCASE)")
        cw.writer(fstream, "#define egtransetCOFLAM egtransetcoflam")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_UNDERSCORE)")
        cw.writer(fstream, "#define egtransetCOFLAM egtransetcoflam_")
        cw.writer(fstream, "#endif")

    # visco coefs
    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
    cw.writer(fstream, "void egtransetCOFLAM(amrex::Real* COFLAM) {")

    for spec in species_info.nonqss_species:
        for i in range(4):
            cw.writer(
                fstream,
                "%s[%d] = %.8E;"
                % ("COFLAM", spec.idx * 4 + i, coflam[spec.idx][3 - i]),
            )

    cw.writer(fstream, "}")


def diffcoefs(fstream, species_info, speciesTransport, do_declarations, NTFit):
    """Write the diffusion coefficients."""
    # REORDERING OF SPECS
    specOrdered = []
    nSpecies = species_info.nSpecies
    for i in range(nSpecies):
        for spec in speciesTransport:
            if spec.idx == i:
                specOrdered.append(spec)
                break

    # compute single constants in g/cm/s
    Na = 6.02214199e23
    # conversion coefs
    AtoCM = 1.0e-8
    DEBYEtoCGS = 1.0e-18
    PATM = 0.1013250000000000e07
    # temperature increment
    dt = (species_info.highT - species_info.lowT) / (NTFit - 1)
    # diff coefs (4 per spec pair)
    cofd = []
    for i, spec1 in enumerate(specOrdered):
        cofd.append([])
        if i != spec1.idx:
            print("Problem in _diffcoefs computation")
            sys.exit(1)
        for j, spec2 in enumerate(specOrdered[0 : i + 1]):
            if j != spec2.idx:
                print("Problem in _diffcoefs computation")
                sys.exit(1)
            # eq. (9)
            sigm = (
                0.5
                * (
                    float(speciesTransport[spec1][2])
                    + float(speciesTransport[spec2][2])
                )
                * AtoCM
            ) * Xi(spec1, spec2, speciesTransport) ** (1.0 / 6.0)
            # eq. (4)
            m_red = (
                spec1.weight
                * spec2.weight
                / (spec1.weight + spec2.weight)
                / Na
            )
            # eq. (8) & (14)
            epsm_k = (
                np.sqrt(
                    float(speciesTransport[spec1][1])
                    * float(speciesTransport[spec2][1])
                )
                * Xi(spec1, spec2, speciesTransport) ** 2.0
            )

            # eq. (15)
            conversion = DEBYEtoCGS * DEBYEtoCGS / cc.kb
            dst = (
                0.5
                * conversion
                * float(speciesTransport[spec1][3])
                * float(speciesTransport[spec2][3])
                / (epsm_k * sigm**3)
            )
            if not Xi_bool(spec1, spec2, speciesTransport):
                dst = 0.0
            # enter the loop on temperature
            spdiffcoef = []
            tlog = []
            for n in range(NTFit):
                t = species_info.lowT + dt * n
                tr = t / epsm_k
                # eq. (3)
                # note: these are "corrected" in CHEMKIN not in CANTERA... we chose not to
                difcoeff = (
                    3.0
                    / 16.0
                    * 1
                    / PATM
                    * (
                        np.sqrt(2.0 * np.pi * t**3 * cc.kb**3 / m_red)
                        / (np.pi * sigm * sigm * om11_CHEMKIN(tr, dst))
                    )
                )
                # log transformation for polyfit
                tlog.append(np.log(t))
                spdiffcoef.append(np.log(difcoeff))

            cofd[i].append(np.polyfit(tlog, spdiffcoef, 3))

    # use the symmetry for upper triangular terms
    # note: starting with this would be preferable (only one bigger loop)
    # note2: or write stuff differently !
    # for i,spec1 in enumerate(specOrdered):
    #    for j,spec2 in enumerate(specOrdered[i+1:]):
    #        cofd[i].append(cofd[spec2.id][spec1.id])

    # header for diffusion coefs
    cw.writer(fstream)
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("Poly fits for the diffusion coefficients, dim NO*KK*KK"),
    )
    if do_declarations:
        cw.writer(fstream, "#if defined(BL_FORT_USE_UPPERCASE)")
        cw.writer(fstream, "#define egtransetCOFD EGTRANSETCOFD")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_LOWERCASE)")
        cw.writer(fstream, "#define egtransetCOFD egtransetcofd")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_UNDERSCORE)")
        cw.writer(fstream, "#define egtransetCOFD egtransetcofd_")
        cw.writer(fstream, "#endif")

    # coefs
    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
    cw.writer(fstream, "void egtransetCOFD(amrex::Real* COFD) {")

    for i, spec1 in enumerate(specOrdered):
        # for j,spec2 in enumerate(specOrdered):
        for j, spec2 in enumerate(specOrdered[0 : i + 1]):
            for k in range(4):
                # cw.writer(fstream,'%s[%d] = %.8E;' % ('COFD', i*nSpecies*4+j*4+k, cofd[j][i][3-k]))
                cw.writer(
                    fstream,
                    "%s[%d] = %.8E;"
                    % (
                        "COFD",
                        i * nSpecies * 4 + j * 4 + k,
                        cofd[i][j][3 - k],
                    ),
                )
        for j, spec2 in enumerate(specOrdered[i + 1 :]):
            for k in range(4):
                cw.writer(
                    fstream,
                    "%s[%d] = %.8E;"
                    % (
                        "COFD",
                        i * nSpecies * 4 + (j + i + 1) * 4 + k,
                        cofd[j + i + 1][i][3 - k],
                    ),
                )

    cw.writer(fstream, "}")


def lightSpecs(fstream, speclist, do_declarations):
    """Write list of specs with small weight, dim NLITE."""
    # header
    cw.writer(fstream)
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("List of specs with small weight, dim NLITE")
    )
    if do_declarations:
        cw.writer(fstream, "#if defined(BL_FORT_USE_UPPERCASE)")
        cw.writer(fstream, "#define egtransetKTDIF EGTRANSETKTDIF")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_LOWERCASE)")
        cw.writer(fstream, "#define egtransetKTDIF egtransetktdif")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_UNDERSCORE)")
        cw.writer(fstream, "#define egtransetKTDIF egtransetktdif_")
        cw.writer(fstream, "#endif")

    # coefs
    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
    cw.writer(fstream, "void egtransetKTDIF(int* KTDIF) {")

    for i in range(len(speclist)):
        cw.writer(fstream, "%s[%d] = %d;" % ("KTDIF", i, speclist[i]))

    cw.writer(fstream, "}")


def thermaldiffratios(
    fstream,
    species_info,
    speciesTransport,
    lightSpecList,
    do_declarations,
    NTFit,
):
    """Write thermal diffusion ratios."""
    nSpecies = species_info.nSpecies
    # This is an overhaul of CHEMKIN version III
    # REORDERING OF SPECS
    specOrdered = []
    for i in range(nSpecies):
        for spec in speciesTransport:
            if spec.idx == i:
                specOrdered.append(spec)
                break

    # compute single constants in g/cm/s
    # conversion coefs
    DEBYEtoCGS = 1.0e-18
    AtoCM = 1.0e-8
    # temperature increment
    dt = (species_info.highT - species_info.lowT) / (NTFit - 1)
    # diff ratios (4 per spec pair involving light species)
    coftd = []
    k = -1
    for i, spec1 in enumerate(specOrdered):
        if i != spec1.idx:
            print("Problem in _thermaldiffratios computation")
            sys.exit(1)
        if spec1.idx in lightSpecList:
            k = k + 1
            if lightSpecList[k] != spec1.idx:
                print("Problem in  _thermaldiffratios computation")
                sys.exit(1)
            coftd.append([])
            epsi = float(speciesTransport[spec1][1]) * cc.kb
            sigi = float(speciesTransport[spec1][2]) * AtoCM
            poli = float(speciesTransport[spec1][4]) * AtoCM * AtoCM * AtoCM
            # eq. (12)
            poliRed = poli / sigi**3
            for j, spec2 in enumerate(specOrdered):
                if j != spec2.idx:
                    print("Problem in _thermaldiffratios computation")
                    sys.exit(1)
                # eq. (53)
                Wji = (spec2.weight - spec1.weight) / (
                    spec1.weight + spec2.weight
                )
                epsj = float(speciesTransport[spec2][1]) * cc.kb
                sigj = float(speciesTransport[spec2][2]) * AtoCM
                dipj = float(speciesTransport[spec2][3]) * DEBYEtoCGS
                # eq. (13)
                dipjRed = dipj / np.sqrt(epsj * sigj**3)
                epsRatio = epsj / epsi
                tse = 1.0 + 0.25 * poliRed * dipjRed**2 * np.sqrt(epsRatio)
                eok = tse**2 * np.sqrt(
                    float(speciesTransport[spec1][1])
                    * float(speciesTransport[spec2][1])
                )
                # enter the loop on temperature
                spthdiffcoef = []
                tTab = []
                for n in range(NTFit):
                    t = species_info.lowT + dt * n
                    tslog = np.log(t) - np.log(eok)
                    # eq. (53)
                    thdifcoeff = (
                        15.0
                        / 2.0
                        * Wji
                        * (2.0 * astar(tslog) + 5.0)
                        * (6.0 * cstar(tslog) - 5.0)
                        / (
                            astar(tslog)
                            * (
                                16.0 * astar(tslog)
                                - 12.0 * bstar(tslog)
                                + 55.0
                            )
                        )
                    )

                    # log transformation for polyfit
                    tTab.append(t)
                    spthdiffcoef.append(thdifcoeff)

                coftd[k].append(np.polyfit(tTab, spthdiffcoef, 3))

    # header for thermal diff ratios
    cw.writer(fstream)
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("Poly fits for thermal diff ratios, dim NO*NLITE*KK"),
    )
    if do_declarations:
        cw.writer(fstream, "#if defined(BL_FORT_USE_UPPERCASE)")
        cw.writer(fstream, "#define egtransetCOFTD EGTRANSETCOFTD")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_LOWERCASE)")
        cw.writer(fstream, "#define egtransetCOFTD egtransetcoftd")
        cw.writer(fstream, "#elif defined(BL_FORT_USE_UNDERSCORE)")
        cw.writer(fstream, "#define egtransetCOFTD egtransetcoftd_")
        cw.writer(fstream, "#endif")

    # visco coefs
    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
    cw.writer(fstream, "void egtransetCOFTD(amrex::Real* COFTD) {")

    for i in range(len(coftd)):
        for j in range(nSpecies):
            for k in range(4):
                cw.writer(
                    fstream,
                    "%s[%d] = %.8E;"
                    % (
                        "COFTD",
                        i * 4 * nSpecies + j * 4 + k,
                        coftd[i][j][3 - k],
                    ),
                )

    cw.writer(fstream, "}")


def astar(tslog):
    """Compute astar."""
    aTab = [
        0.1106910525e01,
        -0.7065517161e-02,
        -0.1671975393e-01,
        0.1188708609e-01,
        0.7569367323e-03,
        -0.1313998345e-02,
        0.1720853282e-03,
    ]

    B = aTab[6]
    for i in range(6):
        B = aTab[5 - i] + B * tslog

    return B


def bstar(tslog):
    """Compute bstar."""
    bTab = [
        0.1199673577e01,
        -0.1140928763e00,
        -0.2147636665e-02,
        0.2512965407e-01,
        -0.3030372973e-02,
        -0.1445009039e-02,
        0.2492954809e-03,
    ]

    B = bTab[6]
    for i in range(6):
        B = bTab[5 - i] + B * tslog

    return B


def cstar(tslog):
    """Compute cstar."""
    cTab = [
        0.8386993788e00,
        0.4748325276e-01,
        0.3250097527e-01,
        -0.1625859588e-01,
        -0.2260153363e-02,
        0.1844922811e-02,
        -0.2115417788e-03,
    ]

    B = cTab[6]
    for i in range(6):
        B = cTab[5 - i] + B * tslog

    return B


def om22_CHEMKIN(tr, dst):
    """Compute OM22."""
    # This is an overhaul of CANTERA version 2.3
    # range of dst
    dstTab = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5]

    # range of tr
    trTab = [
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
        1.2,
        1.4,
        1.6,
        1.8,
        2.0,
        2.5,
        3.0,
        3.5,
        4.0,
        5.0,
        6.0,
        7.0,
        8.0,
        9.0,
        10.0,
        12.0,
        14.0,
        16.0,
        18.0,
        20.0,
        25.0,
        30.0,
        35.0,
        40.0,
        50.0,
        75.0,
        100.0,
    ]

    # tab of omega22 corresp. to (tr, dst)
    # CANTERA
    omegaTab = [
        4.1005,
        4.266,
        4.833,
        5.742,
        6.729,
        8.624,
        10.34,
        11.89,
        3.2626,
        3.305,
        3.516,
        3.914,
        4.433,
        5.57,
        6.637,
        7.618,
        2.8399,
        2.836,
        2.936,
        3.168,
        3.511,
        4.329,
        5.126,
        5.874,
        2.531,
        2.522,
        2.586,
        2.749,
        3.004,
        3.64,
        4.282,
        4.895,
        2.2837,
        2.277,
        2.329,
        2.46,
        2.665,
        3.187,
        3.727,
        4.249,
        2.0838,
        2.081,
        2.13,
        2.243,
        2.417,
        2.862,
        3.329,
        3.786,
        1.922,
        1.924,
        1.97,
        2.072,
        2.225,
        2.614,
        3.028,
        3.435,
        1.7902,
        1.795,
        1.84,
        1.934,
        2.07,
        2.417,
        2.788,
        3.156,
        1.6823,
        1.689,
        1.733,
        1.82,
        1.944,
        2.258,
        2.596,
        2.933,
        1.5929,
        1.601,
        1.644,
        1.725,
        1.838,
        2.124,
        2.435,
        2.746,
        1.4551,
        1.465,
        1.504,
        1.574,
        1.67,
        1.913,
        2.181,
        2.451,
        1.3551,
        1.365,
        1.4,
        1.461,
        1.544,
        1.754,
        1.989,
        2.228,
        1.28,
        1.289,
        1.321,
        1.374,
        1.447,
        1.63,
        1.838,
        2.053,
        1.2219,
        1.231,
        1.259,
        1.306,
        1.37,
        1.532,
        1.718,
        1.912,
        1.1757,
        1.184,
        1.209,
        1.251,
        1.307,
        1.451,
        1.618,
        1.795,
        1.0933,
        1.1,
        1.119,
        1.15,
        1.193,
        1.304,
        1.435,
        1.578,
        1.0388,
        1.044,
        1.059,
        1.083,
        1.117,
        1.204,
        1.31,
        1.428,
        0.99963,
        1.004,
        1.016,
        1.035,
        1.062,
        1.133,
        1.22,
        1.319,
        0.96988,
        0.9732,
        0.983,
        0.9991,
        1.021,
        1.079,
        1.153,
        1.236,
        0.92676,
        0.9291,
        0.936,
        0.9473,
        0.9628,
        1.005,
        1.058,
        1.121,
        0.89616,
        0.8979,
        0.903,
        0.9114,
        0.923,
        0.9545,
        0.9955,
        1.044,
        0.87272,
        0.8741,
        0.878,
        0.8845,
        0.8935,
        0.9181,
        0.9505,
        0.9893,
        0.85379,
        0.8549,
        0.858,
        0.8632,
        0.8703,
        0.8901,
        0.9164,
        0.9482,
        0.83795,
        0.8388,
        0.8414,
        0.8456,
        0.8515,
        0.8678,
        0.8895,
        0.916,
        0.82435,
        0.8251,
        0.8273,
        0.8308,
        0.8356,
        0.8493,
        0.8676,
        0.8901,
        0.80184,
        0.8024,
        0.8039,
        0.8065,
        0.8101,
        0.8201,
        0.8337,
        0.8504,
        0.78363,
        0.784,
        0.7852,
        0.7872,
        0.7899,
        0.7976,
        0.8081,
        0.8212,
        0.76834,
        0.7687,
        0.7696,
        0.7712,
        0.7733,
        0.7794,
        0.7878,
        0.7983,
        0.75518,
        0.7554,
        0.7562,
        0.7575,
        0.7592,
        0.7642,
        0.7711,
        0.7797,
        0.74364,
        0.7438,
        0.7445,
        0.7455,
        0.747,
        0.7512,
        0.7569,
        0.7642,
        0.71982,
        0.72,
        0.7204,
        0.7211,
        0.7221,
        0.725,
        0.7289,
        0.7339,
        0.70097,
        0.7011,
        0.7014,
        0.7019,
        0.7026,
        0.7047,
        0.7076,
        0.7112,
        0.68545,
        0.6855,
        0.6858,
        0.6861,
        0.6867,
        0.6883,
        0.6905,
        0.6932,
        0.67232,
        0.6724,
        0.6726,
        0.6728,
        0.6733,
        0.6743,
        0.6762,
        0.6784,
        0.65099,
        0.651,
        0.6512,
        0.6513,
        0.6516,
        0.6524,
        0.6534,
        0.6546,
        0.61397,
        0.6141,
        0.6143,
        0.6145,
        0.6147,
        0.6148,
        0.6148,
        0.6147,
        0.5887,
        0.5889,
        0.5894,
        0.59,
        0.5903,
        0.5901,
        0.5895,
        0.5885,
    ]

    # First test on tr
    if tr > 75.0:
        omeg12 = (
            0.703
            - 0.146e-2 * tr
            + 0.357e-5 * tr * tr
            - 0.343e-8 * tr * tr * tr
        )
    else:
        # Find tr idx in trTab
        if tr <= 0.2:
            ii = 1
        else:
            ii = 36
        for i in range(1, 37):
            if (tr > trTab[i - 1]) and (tr <= trTab[i]):
                ii = i
                break
        # Find dst idx in dstTab
        if abs(dst) >= 1.0e-5:
            if dst <= 0.25:
                kk = 1
            else:
                kk = 6
            for i in range(1, 7):
                if (dstTab[i - 1] < dst) and (dstTab[i] >= dst):
                    kk = i
                    break
            # Find surrounding values and interpolate
            # First on dst
            vert = np.zeros(3)
            for i in range(3):
                arg = np.zeros(3)
                val = np.zeros(3)
                for k in range(3):
                    arg[k] = dstTab[kk - 1 + k]
                    val[k] = omegaTab[8 * (ii - 1 + i) + (kk - 1 + k)]
                vert[i] = qinterp(dst, arg, val)
            # Second on tr
            arg = np.zeros(3)
            for i in range(3):
                arg[i] = trTab[ii - 1 + i]
            omeg12 = qinterp(tr, arg, vert)
        else:
            arg = np.zeros(3)
            val = np.zeros(3)
            for i in range(3):
                arg[i] = trTab[ii - 1 + i]
                val[i] = omegaTab[8 * (ii - 1 + i)]
            omeg12 = qinterp(tr, arg, val)
    return omeg12


def om11_CHEMKIN(tr, dst):
    """Compute OM11."""
    # This is an overhaul of CANTERA version 2.3
    # range of dst
    dstTab = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5]

    # range of tr
    trTab = [
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
        1.2,
        1.4,
        1.6,
        1.8,
        2.0,
        2.5,
        3.0,
        3.5,
        4.0,
        5.0,
        6.0,
        7.0,
        8.0,
        9.0,
        10.0,
        12.0,
        14.0,
        16.0,
        18.0,
        20.0,
        25.0,
        30.0,
        35.0,
        40.0,
        50.0,
        75.0,
        100.0,
    ]

    # tab of omega11 corresp. to (tr, dst)
    # CANTERA
    omegaTab = [
        4.008,
        4.002,
        4.655,
        5.52,
        6.454,
        8.214,
        9.824,
        11.31,
        3.130,
        3.164,
        3.355,
        3.721,
        4.198,
        5.23,
        6.225,
        7.160,
        2.649,
        2.657,
        2.77,
        3.002,
        3.319,
        4.054,
        4.785,
        5.483,
        2.314,
        2.32,
        2.402,
        2.572,
        2.812,
        3.386,
        3.972,
        4.539,
        2.066,
        2.073,
        2.14,
        2.278,
        2.472,
        2.946,
        3.437,
        3.918,
        1.877,
        1.885,
        1.944,
        2.06,
        2.225,
        2.628,
        3.054,
        3.747,
        1.729,
        1.738,
        1.79,
        1.893,
        2.036,
        2.388,
        2.763,
        3.137,
        1.6122,
        1.622,
        1.67,
        1.76,
        1.886,
        2.198,
        2.535,
        2.872,
        1.517,
        1.527,
        1.572,
        1.653,
        1.765,
        2.044,
        2.35,
        2.657,
        1.44,
        1.45,
        1.49,
        1.564,
        1.665,
        1.917,
        2.196,
        2.4780,
        1.3204,
        1.33,
        1.364,
        1.425,
        1.51,
        1.72,
        1.956,
        2.199,
        1.234,
        1.24,
        1.272,
        1.324,
        1.394,
        1.573,
        1.777,
        1.99,
        1.168,
        1.176,
        1.202,
        1.246,
        1.306,
        1.46,
        1.64,
        1.827,
        1.1166,
        1.124,
        1.146,
        1.185,
        1.237,
        1.372,
        1.53,
        1.7,
        1.075,
        1.082,
        1.102,
        1.135,
        1.181,
        1.3,
        1.441,
        1.592,
        1.0006,
        1.005,
        1.02,
        1.046,
        1.08,
        1.17,
        1.278,
        1.397,
        0.95,
        0.9538,
        0.9656,
        0.9852,
        1.012,
        1.082,
        1.168,
        1.265,
        0.9131,
        0.9162,
        0.9256,
        0.9413,
        0.9626,
        1.019,
        1.09,
        1.17,
        0.8845,
        0.8871,
        0.8948,
        0.9076,
        0.9252,
        0.972,
        1.03,
        1.098,
        0.8428,
        0.8446,
        0.850,
        0.859,
        0.8716,
        0.9053,
        0.9483,
        0.9984,
        0.813,
        0.8142,
        0.8183,
        0.825,
        0.8344,
        0.8598,
        0.8927,
        0.9316,
        0.7898,
        0.791,
        0.794,
        0.7993,
        0.8066,
        0.8265,
        0.8526,
        0.8836,
        0.7711,
        0.772,
        0.7745,
        0.7788,
        0.7846,
        0.8007,
        0.822,
        0.8474,
        0.7555,
        0.7562,
        0.7584,
        0.7619,
        0.7667,
        0.78,
        0.7976,
        0.8189,
        0.7422,
        0.743,
        0.7446,
        0.7475,
        0.7515,
        0.7627,
        0.7776,
        0.796,
        0.72022,
        0.7206,
        0.722,
        0.7241,
        0.7271,
        0.7354,
        0.7464,
        0.76,
        0.7025,
        0.703,
        0.704,
        0.7055,
        0.7078,
        0.7142,
        0.7228,
        0.7334,
        0.68776,
        0.688,
        0.6888,
        0.6901,
        0.6919,
        0.697,
        0.704,
        0.7125,
        0.6751,
        0.6753,
        0.676,
        0.677,
        0.6785,
        0.6827,
        0.6884,
        0.6955,
        0.664,
        0.6642,
        0.6648,
        0.6657,
        0.6669,
        0.6704,
        0.6752,
        0.681,
        0.6414,
        0.6415,
        0.6418,
        0.6425,
        0.6433,
        0.6457,
        0.649,
        0.653,
        0.6235,
        0.6236,
        0.6239,
        0.6243,
        0.6249,
        0.6267,
        0.629,
        0.632,
        0.60882,
        0.6089,
        0.6091,
        0.6094,
        0.61,
        0.6112,
        0.613,
        0.6154,
        0.5964,
        0.5964,
        0.5966,
        0.597,
        0.5972,
        0.5983,
        0.600,
        0.6017,
        0.5763,
        0.5763,
        0.5764,
        0.5766,
        0.5768,
        0.5775,
        0.5785,
        0.58,
        0.5415,
        0.5415,
        0.5416,
        0.5416,
        0.5418,
        0.542,
        0.5424,
        0.543,
        0.518,
        0.518,
        0.5182,
        0.5184,
        0.5184,
        0.5185,
        0.5186,
        0.5187,
    ]

    # First test on tr
    if tr > 75.0:
        omeg12 = (
            0.623
            - 0.136e-2 * tr
            + 0.346e-5 * tr * tr
            - 0.343e-8 * tr * tr * tr
        )
    else:
        # Find tr idx in trTab
        if tr <= 0.2:
            ii = 1
        else:
            ii = 36
        for i in range(1, 37):
            if (tr > trTab[i - 1]) and (tr <= trTab[i]):
                ii = i
                break
        # Find dst idx in dstTab
        if abs(dst) >= 1.0e-5:
            if dst <= 0.25:
                kk = 1
            else:
                kk = 6
            for i in range(1, 7):
                if (dstTab[i - 1] < dst) and (dstTab[i] >= dst):
                    kk = i
                    break
            # Find surrounding values and interpolate
            # First on dst
            vert = np.zeros(3)
            for i in range(3):
                arg = np.zeros(3)
                val = np.zeros(3)
                for k in range(3):
                    arg[k] = dstTab[kk - 1 + k]
                    val[k] = omegaTab[8 * (ii - 1 + i) + (kk - 1 + k)]
                vert[i] = qinterp(dst, arg, val)
            # Second on tr
            arg = np.zeros(3)
            for i in range(3):
                arg[i] = trTab[ii - 1 + i]
            omeg12 = qinterp(tr, arg, vert)
        else:
            arg = np.zeros(3)
            val = np.zeros(3)
            for i in range(3):
                arg[i] = trTab[ii - 1 + i]
                val[i] = omegaTab[8 * (ii - 1 + i)]
            omeg12 = qinterp(tr, arg, val)
    return omeg12


def qinterp(x0, x, y):
    """Compute interpolated value."""
    val1 = y[0] + (x0 - x[0]) * (y[1] - y[0]) / (x[1] - x[0])
    val2 = y[1] + (x0 - x[1]) * (y[2] - y[1]) / (x[2] - x[1])
    fac1 = (x0 - x[0]) / (x[1] - x[0]) / 2.0
    fac2 = (x[2] - x0) / (x[2] - x[1]) / 2.0
    if x0 >= x[1]:
        val = (val1 * fac2 + val2) / (1.0 + fac2)
    else:
        val = (val1 + val2 * fac1) / (1.0 + fac1)
    return val


def getCVdRspecies(mechanism, t, species):
    """Get the parameters of a thermo model."""
    model = mechanism.species(species.name).thermo
    if not model.n_coeffs == 15:
        print("Unsupported thermo model.")
        sys.exit(1)

    mid = model.coeffs[0]
    highRange = model.coeffs[1:8]
    lowRange = model.coeffs[8:15]

    if t < mid:
        parameters = lowRange
    else:
        parameters = highRange

    return (
        (parameters[0] - 1.0)
        + parameters[1] * t
        + parameters[2] * t * t
        + parameters[3] * t * t * t
        + parameters[4] * t * t * t * t
    )


def Fcorr(t, eps_k):
    """Compute Fcorr value."""
    thtwo = 3.0 / 2.0
    return (
        1
        + np.pi ** (thtwo) / 2.0 * np.sqrt((eps_k / t))
        + (np.pi**2 / 4.0 + 2.0) * ((eps_k / t))
        + ((np.pi * eps_k / t)) ** (thtwo)
    )


def Xi(spec1, spec2, speciesTransport):
    """Compute Xi."""
    dipmin = 1e-20
    # 1 is polar, 2 is nonpolar
    # err in eq. (11) ?
    if (float(speciesTransport[spec2][3]) < dipmin) and (
        float(speciesTransport[spec1][3]) > dipmin
    ):
        xi = 1.0 + 1.0 / 4.0 * redPol(spec2, speciesTransport) * redDip(
            spec1, speciesTransport
        ) * redDip(spec1, speciesTransport) * np.sqrt(
            float(speciesTransport[spec1][1])
            / float(speciesTransport[spec2][1])
        )
    # 1 is nonpolar, 2 is polar
    elif (float(speciesTransport[spec2][3]) > dipmin) and (
        float(speciesTransport[spec1][3]) < dipmin
    ):
        xi = 1.0 + 1.0 / 4.0 * redPol(spec1, speciesTransport) * redDip(
            spec2, speciesTransport
        ) * redDip(spec2, speciesTransport) * np.sqrt(
            float(speciesTransport[spec2][1])
            / float(speciesTransport[spec1][1])
        )
    # normal case, either both polar or both nonpolar
    else:
        xi = 1.0

    return xi


def Xi_bool(spec1, spec2, speciesTransport):
    """Compute the boolean of Xi."""
    dipmin = 1e-20
    # 1 is polar, 2 is nonpolar
    # err in eq. (11) ?
    if (float(speciesTransport[spec2][3]) < dipmin) and (
        float(speciesTransport[spec1][3]) > dipmin
    ):
        xi_b = False
    # 1 is nonpolar, 2 is polar
    elif (float(speciesTransport[spec2][3]) > dipmin) and (
        float(speciesTransport[spec1][3]) < dipmin
    ):
        xi_b = False
    # normal case, either both polar or both nonpolar
    else:
        xi_b = True

    return xi_b


def redPol(spec, speciesTransport):
    """Compute polarization value."""
    return (
        float(speciesTransport[spec][4])
        / float(speciesTransport[spec][2]) ** 3.0
    )


def redDip(spec, speciesTransport):
    """Compute dipole value."""
    # compute single constants in g/cm/s
    # conversion coefs
    AtoCM = 1.0e-8
    DEBYEtoCGS = 1.0e-18
    convert = DEBYEtoCGS / np.sqrt(cc.kb * AtoCM**3.0)
    return (
        convert
        * float(speciesTransport[spec][3])
        / np.sqrt(
            float(speciesTransport[spec][1])
            * float(speciesTransport[spec][2]) ** 3.0
        )
    )


def getCriticalParameters(fstream, mechanism, species_info):
    """Write the critical parameters."""
    TabulatedCriticalParams = {
        "H2": {
            "Tci": 33.145,
            "Pci": 12.964,
            "wt": 2.01588,
            "acentric_factor": -0.219,
        },
        "O2": {
            "Tci": 154.581,
            "Pci": 50.4304658,
            "wt": 31.9988,
            "acentric_factor": 0.0222,
        },
        "H2O": {
            "Tci": 647.096,
            "Pci": 220.640,
            "wt": 18.015340,
            "acentric_factor": 0.3443,
        },
        "N2": {
            "Tci": 126.192,
            "Pci": 33.958,
            "wt": 28.013400,
            "acentric_factor": 0.0372,
        },
        "CH4": {
            "Tci": 190.56,
            "Pci": 45.99,
            "wt": 16.043030,
            "acentric_factor": 0.011,
        },
        "C2H6": {
            "Tci": 305.32,
            "Pci": 48.72,
            "wt": 30.070120,
            "acentric_factor": 0.099,
        },
        "C3H8": {
            "Tci": 369.83,
            "Pci": 42.48,
            "wt": 44.097210,
            "acentric_factor": 0.152,
        },
        "CO2": {
            "Tci": 304.12,
            "Pci": 73.74,
            "wt": 44.009950,
            "acentric_factor": 0.225,
        },
        "He": {
            "Tci": 5.1953,
            "Pci": 2.2746,
            "wt": 4.002602,
            "acentric_factor": -0.382,
        },
        "CO": {
            "Tci": 132.85,
            "Pci": 34.94,
            "wt": 28.010,
            "acentric_factor": 0.045,
        },
        "AR": {
            "Tci": 150.86,
            "Pci": 48.98,
            "wt": 39.948,
            "acentric_factor": -0.002,
        },
        "NO": {
            "Tci": 180.0,
            "Pci": 64.80,
            "wt": 30.006,
            "acentric_factor": 0.582,
        },
        "CH3OH": {
            "Tci": 512.64,
            "Pci": 80.97,
            "wt": 32.042,
            "acentric_factor": 0.565,
        },
        "C2H2": {
            "Tci": 308.30,
            "Pci": 61.14,
            "wt": 26.038,
            "acentric_factor": 0.189,
        },
        "C2H4": {
            "Tci": 282.34,
            "Pci": 50.41,
            "wt": 28.054,
            "acentric_factor": 0.087,
        },
        "N2O": {
            "Tci": 309.60,
            "Pci": 72.55,
            "wt": 44.013,
            "acentric_factor": 0.162,
        },
    }

    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("compute the critical parameters for each species")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void"
        " GET_CRITPARAMS(amrex::Real *  Tci, amrex::Real *  ai, amrex::Real * "
        " bi, amrex::Real *  acentric_i)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream)

    cw.writer(fstream, "amrex::Real   EPS[%d];" % nSpecies)
    cw.writer(fstream, "amrex::Real   SIG[%d];" % nSpecies)
    cw.writer(fstream, "amrex::Real    wt[%d];" % nSpecies)
    cw.writer(fstream, "amrex::Real avogadro = 6.02214199e23;")
    cw.writer(
        fstream, "amrex::Real boltzmann = 1.3806503e-16; //we work in CGS"
    )
    cw.writer(fstream, "amrex::Real Rcst = 83.144598; //in bar [CGS] !")

    cw.writer(fstream)

    cw.writer(fstream, "egtransetEPS(EPS);")
    cw.writer(fstream, "egtransetSIG(SIG);")
    cw.writer(fstream, "get_mw(wt);")

    for species in species_info.nonqss_species:
        if species.name in TabulatedCriticalParams:
            cw.writer(fstream)
            cw.writer(
                fstream,
                cw.comment("species %d: %s" % (species.idx, species.name)),
            )
            cw.writer(fstream, cw.comment("Imported from NIST"))
            cw.writer(
                fstream,
                "Tci[%d] = %f ; "
                % (
                    species.idx,
                    TabulatedCriticalParams[species.name]["Tci"],
                ),
            )
            cw.writer(
                fstream,
                "ai[%d] = 1e6 * 0.42748 * Rcst * Rcst * Tci[%d] * Tci[%d] /"
                " (%f * %f * %f); "
                % (
                    species.idx,
                    species.idx,
                    species.idx,
                    TabulatedCriticalParams[species.name]["wt"],
                    TabulatedCriticalParams[species.name]["wt"],
                    TabulatedCriticalParams[species.name]["Pci"],
                ),
            )
            cw.writer(
                fstream,
                "bi[%d] = 0.08664 * Rcst * Tci[%d] / (%f * %f); "
                % (
                    species.idx,
                    species.idx,
                    TabulatedCriticalParams[species.name]["wt"],
                    TabulatedCriticalParams[species.name]["Pci"],
                ),
            )
            cw.writer(
                fstream,
                "acentric_i[%d] = %f ;"
                % (
                    species.idx,
                    TabulatedCriticalParams[species.name]["acentric_factor"],
                ),
            )
        else:

            cw.writer(fstream)
            cw.writer(
                fstream,
                cw.comment("species %d: %s" % (species.idx, species.name)),
            )
            cw.writer(
                fstream,
                "Tci[%d] = 1.316 * EPS[%d] ; " % (species.idx, species.idx),
            )
            cw.writer(
                fstream,
                "ai[%d] = (5.55 * avogadro * avogadro * EPS[%d]*boltzmann *"
                " pow(1e-8*SIG[%d],3.0) ) / (wt[%d] * wt[%d]); "
                % (
                    species.idx,
                    species.idx,
                    species.idx,
                    species.idx,
                    species.idx,
                ),
            )
            cw.writer(
                fstream,
                "bi[%d] = 0.855 * avogadro * pow(1e-8*SIG[%d],3.0) /"
                " (wt[%d]); " % (species.idx, species.idx, species.idx),
            )
            cw.writer(fstream, "acentric_i[%d] = 0.0 ;" % (species.idx))

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")