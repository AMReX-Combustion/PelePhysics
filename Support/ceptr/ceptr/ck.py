import ceptr.constants as cc
import ceptr.writer as cw


def ckawt(fstream, mechanism):
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get atomic weight for all elements"))
    cw.writer(fstream, "void CKAWT" + cc.sym + "( amrex::Real *  awt)")
    cw.writer(fstream, "{")
    cw.writer(fstream, "atomicWeight(awt);")
    cw.writer(fstream, "}")


def ckncf(fstream, mechanism, species_info):
    nElement = mechanism.n_elements
    nSpecies = species_info.nSpecies

    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns the elemental composition "))
    cw.writer(fstream, cw.comment("of the speciesi (mdim is num of elements)"))
    cw.writer(fstream, "void CKNCF" + cc.sym + "(int * ncf)")
    cw.writer(fstream, "{")
    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(fstream, "int kd = %d; " % (nElement))
    cw.writer(fstream, cw.comment("Zero ncf"))
    cw.writer(fstream, "for (id = 0; id < kd * %d; ++ id) {" % (nSpecies))
    cw.writer(fstream, " ncf[id] = 0; ")
    cw.writer(fstream, "}")

    cw.writer(fstream)
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(fstream, cw.comment("%s" % species.name))
        for elem, coef in mechanism.species(sp).composition.items():
            cw.writer(
                fstream,
                "ncf[ %d * kd + %d ] = %d; "
                % (
                    species_info.ordered_idx_map[sp],
                    mechanism.element_index(elem),
                    coef,
                )
                + cw.comment("%s" % elem),
            )
        cw.writer(fstream)
    cw.writer(fstream, "}")


def cksyme_str(fstream, mechanism, species_info):
    nElement = mechanism.n_elements
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the vector of strings of element names")
    )
    cw.writer(
        fstream,
        "void CKSYME_STR" + cc.sym + "(amrex::Vector<std::string>& ename)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "ename.resize(%d);" % nElement)
    for elem in mechanism.element_names:
        cw.writer(
            fstream,
            'ename[%d] = "%s";' % (mechanism.element_index(elem), elem),
        )
    cw.writer(fstream, "}")


def cksyms_str(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the vector of strings of species names")
    )
    cw.writer(
        fstream,
        "void CKSYMS_STR" + cc.sym + "(amrex::Vector<std::string>& kname)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "kname.resize(%d);" % nSpecies)
    for species in species_info.nonqss_species_list:
        cw.writer(
            fstream,
            'kname[%d] = "%s";'
            % (species_info.ordered_idx_map[species], species),
        )
    cw.writer(fstream, "}")


def ckindx(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("A few mechanism parameters"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKINDX"
        + cc.sym
        + "(int * mm, int * kk, int * ii, int * nfit)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "*mm = %d;" % mechanism.n_elements)
    cw.writer(fstream, "*kk = %d;" % mechanism.n_species)
    cw.writer(fstream, "*ii = %d;" % mechanism.n_reactions)
    cw.writer(
        fstream, "*nfit = -1; " + cw.comment("Why do you need this anyway ? ")
    )
    cw.writer(fstream, "}")


def ckrp(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(fstream, cw.comment(" Returns R, Rc, Patm"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRP"
        + cc.sym
        + "(amrex::Real *  ru, amrex::Real *  ruc, amrex::Real *  pa)",
    )
    cw.writer(fstream, "{")
    cw.writer(
        fstream,
        " *ru  = %1.14e; "
        % ((cc.R * cc.ureg.mole * cc.ureg.kelvin / cc.ureg.erg)).m,
    )
    cw.writer(
        fstream,
        " *ruc = %.20f; "
        % (cc.Rc * (cc.ureg.mole * cc.ureg.kelvin / cc.ureg.cal)),
    )
    cw.writer(fstream, " *pa  = %g; " % (cc.Patm))
    cw.writer(fstream, "}")


def ckcpbl(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the mean specific heat at CP (Eq. 33)")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPBL"
        + cc.sym
        + "(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  cpbl)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(fstream, "amrex::Real result = 0; ")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real cpor[%d]; " % nSpecies + cw.comment(" temporary storage"),
    )

    # call routine
    cw.writer(fstream, "cp_R(cpor, tc);")

    # dot product
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("perform dot product"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)

    cw.writer(fstream, "result += x[id]*cpor[id];")

    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(
        fstream,
        "*cpbl = result * %1.14e;"
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m,
    )
    cw.writer(fstream, "}")


def ckcpbs(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the mean specific heat at CP (Eq. 34)")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPBS"
        + cc.sym
        + "(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  cpbs)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "amrex::Real result = 0; ")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real cpor[%d], tresult[%d]; " % (nSpecies, nSpecies)
        + cw.comment(" temporary storage"),
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # call routine
    cw.writer(fstream, "cp_R(cpor, tc);")
    cw.writer(fstream)

    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")

    cw.writer(fstream, "tresult[i] = cpor[i]*y[i]*imw[i];")

    cw.writer(fstream, "")
    cw.writer(fstream, "}")
    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")

    cw.writer(fstream, "result += tresult[i];")

    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(
        fstream,
        "*cpbs = result * %1.14e;"
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m,
    )
    cw.writer(fstream, "}")


def ckcvbl(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the mean specific heat at CV (Eq. 35)")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBL"
        + cc.sym
        + "(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  cvbl)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(fstream, "amrex::Real result = 0; ")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real cvor[%d]; " % nSpecies + cw.comment(" temporary storage"),
    )

    # call routine
    cw.writer(fstream, "cv_R(cvor, tc);")

    # dot product
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("perform dot product"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)

    cw.writer(fstream, "result += x[id]*cvor[id];")

    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(
        fstream,
        "*cvbl = result * %1.14e;"
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m,
    )
    cw.writer(fstream, "}")


def ckcvbs(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the mean specific heat at CV (Eq. 36)")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBS"
        + cc.sym
        + "(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  cvbs)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "amrex::Real result = 0; ")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real cvor[%d]; " % nSpecies + cw.comment(" temporary storage"),
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # call routine
    cw.writer(fstream, "cv_R(cvor, tc);")
    cw.writer(fstream)

    # do dot product
    cw.writer(fstream, cw.comment("multiply by y/molecularweight"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "result += cvor[%d]*y[%d]*imw[%d]; "
            % (spec_idx, spec_idx, spec_idx)
            + cw.comment("%s" % species.name),
        )

    cw.writer(fstream)
    cw.writer(
        fstream,
        "*cvbs = result * %1.14e;"
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m,
    )
    cw.writer(fstream, "}")


def ckhbml(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("Returns the mean enthalpy of the mixture in molar units"),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBML"
        + cc.sym
        + "(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  hbml)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(fstream, "amrex::Real result = 0; ")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real hml[%d]; " % nSpecies + cw.comment(" temporary storage"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("R*T"),
    )

    # call routine
    cw.writer(fstream, "speciesEnthalpy(hml, tc);")

    # dot product
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("perform dot product"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)

    cw.writer(fstream, "result += x[id]*hml[id];")

    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "*hbml = result * RT;")
    cw.writer(fstream, "}")


def ckhbms(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns mean enthalpy of mixture in mass units")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBMS"
        + cc.sym
        + "(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  hbms)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "amrex::Real result = 0;")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real hml[%d], tmp[%d]; " % (nSpecies, nSpecies)
        + cw.comment(" temporary storage"),
    )

    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("R*T"),
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # call routine
    cw.writer(fstream, "speciesEnthalpy(hml, tc);")
    cw.writer(fstream)

    cw.writer(fstream, "int id;")
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)

    cw.writer(fstream, "tmp[id] = y[id]*hml[id]*imw[id];")

    cw.writer(fstream, "}")
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)

    cw.writer(fstream, "result += tmp[id];")

    cw.writer(fstream, "}")

    cw.writer(fstream)
    # finally, multiply by RT
    cw.writer(fstream, "*hbms = result * RT;")
    cw.writer(fstream, "}")


def ckubml(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get mean internal energy in molar units"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUBML"
        + cc.sym
        + "(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  ubml)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(fstream, "amrex::Real result = 0; ")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real uml[%d]; " % nSpecies
        + cw.comment(" temporary energy array"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("R*T"),
    )

    # call routine
    cw.writer(fstream, "speciesInternalEnergy(uml, tc);")

    # dot product
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("perform dot product"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "result += x[id]*uml[id];")
    cw.writer(fstream, "}")
    cw.writer(fstream)
    cw.writer(fstream, "*ubml = result * RT;")
    cw.writer(fstream, "}")


def ckubms(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get mean internal energy in mass units"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUBMS"
        + cc.sym
        + "(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  ubms)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "amrex::Real result = 0;")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real ums[%d]; " % nSpecies
        + cw.comment(" temporary energy array"),
    )

    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("R*T"),
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # call routine
    cw.writer(fstream, "speciesInternalEnergy(ums, tc);")
    cw.writer(fstream)

    # convert e/RT to e with mass units
    cw.writer(fstream, cw.comment("perform dot product + scaling by wt"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "result += y[%d]*ums[%d]*imw[%d]; "
            % (spec_idx, spec_idx, spec_idx)
            + cw.comment("%s" % species.name),
        )

    cw.writer(fstream)
    # finally, multiply by RT
    cw.writer(fstream, "*ubms = result * RT;")
    cw.writer(fstream, "}")


def cksbml(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get mixture entropy in molar units"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBML"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real"
        " *  sbml)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(fstream, "amrex::Real result = 0; ")

    # get temperature cache
    cw.writer(
        fstream,
        cw.comment(
            "Log of normalized pressure in cgs units dynes/cm^2 by Patm"
        ),
    )
    cw.writer(fstream, "amrex::Real logPratio = log ( *P / 1013250.0 ); ")
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real sor[%d]; " % nSpecies + cw.comment(" temporary storage"),
    )

    # call routine
    cw.writer(fstream, "speciesEntropy(sor, tc);")

    # Equation 42
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute Eq 42"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)

    cw.writer(
        fstream,
        "result += x[id]*(sor[id]-log((x[id]+%g))-logPratio);" % cc.smallnum,
    )

    cw.writer(fstream, "}")

    cw.writer(fstream)

    cw.writer(
        fstream,
        "*sbml = result * %1.14e;"
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m,
    )
    cw.writer(fstream, "}")


def cksbms(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get mixture entropy in mass units"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBMS"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real"
        " *  sbms)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "amrex::Real result = 0; ")

    # get temperature cache
    cw.writer(
        fstream,
        cw.comment(
            "Log of normalized pressure in cgs units dynes/cm^2 by Patm"
        ),
    )
    cw.writer(fstream, "amrex::Real logPratio = log ( *P / 1013250.0 ); ")
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real sor[%d]; " % nSpecies + cw.comment(" temporary storage"),
    )
    cw.writer(
        fstream,
        "amrex::Real x[%d]; " % nSpecies
        + cw.comment(" need a ytx conversion"),
    )

    cw.writer(
        fstream,
        "amrex::Real YOW = 0; " + cw.comment("See Eq 4, 6 in CK Manual"),
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # compute inverse of mean molecular weight first (eq 3)
    cw.writer(
        fstream, cw.comment("Compute inverse of mean molecular wt first")
    )
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "YOW += y[%d]*imw[%d]; " % (spec_idx, spec_idx)
            + cw.comment("%s" % species.name),
        )

    # now to ytx
    cw.writer(fstream, cw.comment("Now compute y to x conversion"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "x[%d] = y[%d]/(%f*YOW); " % (spec_idx, spec_idx, species.weight),
        )

    # call routine
    cw.writer(fstream, "speciesEntropy(sor, tc);")

    # Equation 42 and 43
    cw.writer(fstream, cw.comment("Perform computation in Eq 42 and 43"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        cw.writer(
            fstream,
            "result += x[%d]*(sor[%d]-log((x[%d]+%g))-logPratio);"
            % (spec_idx, spec_idx, spec_idx, cc.smallnum),
        )

    cw.writer(fstream, cw.comment("Scale by R/W"))
    cw.writer(
        fstream,
        "*sbms = result * %1.14e * YOW;"
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m,
    )
    cw.writer(fstream, "}")


def ckgbml(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns mean gibbs free energy in molar units")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKGBML"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real"
        " *  gbml)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(fstream, "amrex::Real result = 0; ")

    # get temperature cache
    cw.writer(
        fstream,
        cw.comment(
            "Log of normalized pressure in cgs units dynes/cm^2 by Patm"
        ),
    )
    cw.writer(fstream, "amrex::Real logPratio = log ( *P / 1013250.0 ); ")
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("R*T"),
    )
    cw.writer(
        fstream,
        "amrex::Real gort[%d]; " % nSpecies + cw.comment(" temporary storage"),
    )

    # call routine
    cw.writer(fstream, cw.comment("Compute g/RT"))
    cw.writer(fstream, "gibbs(gort, tc);")

    # Equation 44
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute Eq 44"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(
        fstream,
        "result += x[id]*(gort[id]+log((x[id]+%g))+logPratio);" % cc.smallnum,
    )
    cw.writer(fstream, "}")
    cw.writer(fstream)
    cw.writer(fstream, "*gbml = result * RT;")
    cw.writer(fstream, "}")


def ckgbms(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns mixture gibbs free energy in mass units")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKGBMS"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real"
        " *  gbms)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "amrex::Real result = 0; ")

    # get temperature cache
    cw.writer(
        fstream,
        cw.comment(
            "Log of normalized pressure in cgs units dynes/cm^2 by Patm"
        ),
    )
    cw.writer(fstream, "amrex::Real logPratio = log ( *P / 1013250.0 ); ")
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("R*T"),
    )
    cw.writer(
        fstream,
        "amrex::Real gort[%d]; " % nSpecies + cw.comment(" temporary storage"),
    )
    cw.writer(
        fstream,
        "amrex::Real x[%d]; " % nSpecies
        + cw.comment(" need a ytx conversion"),
    )

    cw.writer(
        fstream,
        "amrex::Real YOW = 0; " + cw.comment("To hold 1/molecularweight"),
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # compute inverse of mean molecular weight first (eq 3)
    cw.writer(
        fstream, cw.comment("Compute inverse of mean molecular wt first")
    )
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "YOW += y[%d]*imw[%d]; " % (spec_idx, spec_idx)
            + cw.comment("%s" % species.name),
        )

    # now to ytx
    cw.writer(fstream, cw.comment("Now compute y to x conversion"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "x[%d] = y[%d]/(%f*YOW); " % (spec_idx, spec_idx, species.weight),
        )

    # call routine
    cw.writer(fstream, "gibbs(gort, tc);")

    # Equation 42 and 43
    cw.writer(fstream, cw.comment("Perform computation in Eq 44"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        cw.writer(
            fstream,
            "result += x[%d]*(gort[%d]+log((x[%d]+%g))+logPratio);"
            % (spec_idx, spec_idx, spec_idx, cc.smallnum),
        )

    cw.writer(fstream, cw.comment("Scale by RT/W"))
    cw.writer(fstream, "*gbms = result * RT * YOW;")
    cw.writer(fstream, "}")


def ckabml(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("Returns mean helmholtz free energy in molar units"),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKABML"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real"
        " *  abml)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(fstream, "amrex::Real result = 0; ")

    # get temperature cache
    cw.writer(
        fstream,
        cw.comment(
            "Log of normalized pressure in cgs units dynes/cm^2 by Patm"
        ),
    )
    cw.writer(fstream, "amrex::Real logPratio = log ( *P / 1013250.0 ); ")
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("R*T"),
    )
    cw.writer(
        fstream,
        "amrex::Real aort[%d]; " % nSpecies + cw.comment(" temporary storage"),
    )

    # call routine
    cw.writer(fstream, cw.comment("Compute g/RT"))
    cw.writer(fstream, "helmholtz(aort, tc);")

    # Equation 44
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute Eq 44"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)

    cw.writer(
        fstream,
        "result += x[id]*(aort[id]+log((x[id]+%g))+logPratio);" % cc.smallnum,
    )

    cw.writer(fstream, "}")
    cw.writer(fstream)
    cw.writer(fstream, "*abml = result * RT;")
    cw.writer(fstream, "}")


def ckabms(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("Returns mixture helmholtz free energy in mass units"),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKABMS"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real"
        " *  abms)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "amrex::Real result = 0; ")

    # get temperature cache
    cw.writer(
        fstream,
        cw.comment(
            "Log of normalized pressure in cgs units dynes/cm^2 by Patm"
        ),
    )
    cw.writer(fstream, "amrex::Real logPratio = log ( *P / 1013250.0 ); ")
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("R*T"),
    )
    cw.writer(
        fstream,
        "amrex::Real aort[%d]; " % nSpecies + cw.comment(" temporary storage"),
    )
    cw.writer(
        fstream,
        "amrex::Real x[%d]; " % nSpecies
        + cw.comment(" need a ytx conversion"),
    )

    cw.writer(
        fstream,
        "amrex::Real YOW = 0; " + cw.comment("To hold 1/molecularweight"),
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # compute inverse of mean molecular weight first (eq 3)
    cw.writer(
        fstream, cw.comment("Compute inverse of mean molecular wt first")
    )
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "YOW += y[%d]*imw[%d]; " % (spec_idx, spec_idx)
            + cw.comment("%s" % species.name),
        )

    # now to ytx
    cw.writer(fstream, cw.comment("Now compute y to x conversion"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "x[%d] = y[%d]/(%f*YOW); " % (spec_idx, spec_idx, species.weight),
        )

    # call routine
    cw.writer(fstream, "helmholtz(aort, tc);")

    # Equation 42 and 43
    cw.writer(fstream, cw.comment("Perform computation in Eq 44"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        cw.writer(
            fstream,
            "result += x[%d]*(aort[%d]+log((x[%d]+%g))+logPratio);"
            % (spec_idx, spec_idx, spec_idx, cc.smallnum),
        )

    cw.writer(fstream, cw.comment("Scale by RT/W"))
    cw.writer(fstream, "*abms = result * RT * YOW;")

    cw.writer(fstream, "}")


def ckpx(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute P = rhoRT/W(x)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPX"
        + cc.sym
        + "(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x,"
        " amrex::Real *  P)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "amrex::Real XW = 0;" + cw.comment(" To hold mean molecular wt"),
    )

    # molecular weights of all species
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "XW += x[%d]*%f; " % (spec_idx, species.weight)
            + cw.comment("%s" % species.name),
        )

    cw.writer(
        fstream,
        "*P = *rho * %1.14e * (*T) / XW; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("P = rho*R*T/W"),
    )

    cw.writer(fstream)
    cw.writer(fstream, "return;")

    cw.writer(fstream, "}")


def ckpy(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute P = rhoRT/W(y)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPY"
        + cc.sym
        + "(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, "
        " amrex::Real *  P)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream, "amrex::Real YOW = 0;" + cw.comment(" for computing mean MW")
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # molecular weights of all species
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "YOW += y[%d]*imw[%d]; " % (spec_idx, spec_idx)
            + cw.comment("%s" % species.name),
        )

    cw.comment("YOW holds the reciprocal of the mean molecular wt")
    cw.writer(
        fstream,
        "*P = *rho * %1.14e * (*T) * YOW; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("P = rho*R*T/W"),
    )

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckpc(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute P = rhoRT/W(c)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPC"
        + cc.sym
        + "(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  c, "
        " amrex::Real *  P)",
    )

    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(fstream, cw.comment("See Eq 5 in CK Manual"))
    cw.writer(fstream, "amrex::Real W = 0;")
    cw.writer(fstream, "amrex::Real sumC = 0;")

    # molecular weights of all species
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "W += c[%d]*%f; " % (spec_idx, species.weight)
            + cw.comment("%s" % species.name),
        )

    cw.writer(fstream)
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "sumC += c[id];")
    cw.writer(fstream, "}")

    cw.comment("W/sumC holds the mean molecular wt")
    cw.writer(
        fstream,
        "*P = *rho * %1.14e * (*T) * sumC / W; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("P = rho*R*T/W"),
    )

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckrhox(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute rho = PW(x)/RT"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOX"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real"
        " *  rho)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "amrex::Real XW = 0;" + cw.comment(" To hold mean molecular wt"),
    )

    # molecular weights of all species
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "XW += x[%d]*%f; " % (spec_idx, species.weight)
            + cw.comment("%s" % species.name),
        )

    cw.writer(
        fstream,
        "*rho = *P * XW / (%1.14e * (*T)); "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("rho = P*W/(R*T)"),
    )

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckrhoy(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute rho = P*W(y)/RT"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOY"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real"
        " *  rho)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "amrex::Real YOW = 0;")
    cw.writer(fstream, "amrex::Real tmp[%d];" % (nSpecies))
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)
    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "tmp[i] = y[i]*imw[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "YOW += tmp[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "")
    cw.writer(
        fstream,
        "*rho = *P / (%1.14e * (*T) * YOW);"
        % ((cc.R * cc.ureg.mole * cc.ureg.kelvin / cc.ureg.erg)).m
        + cw.comment("rho = P*W/(R*T)"),
    )
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckrhoc(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute rho = P*W(c)/(R*T)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOC"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  c,  amrex::Real"
        " *  rho)",
    )

    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(fstream, cw.comment("See Eq 5 in CK Manual"))
    cw.writer(fstream, "amrex::Real W = 0;")
    cw.writer(fstream, "amrex::Real sumC = 0;")

    # molecular weights of all species
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "W += c[%d]*%f; " % (spec_idx, species.weight)
            + cw.comment("%s" % species.name),
        )

    cw.writer(fstream)
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "sumC += c[id];")
    cw.writer(fstream, "}")

    cw.comment("W/sumC holds the mean molecular wt")
    cw.writer(
        fstream,
        "*rho = *P * W / (sumC * (*T) * %1.14e); "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m)
        + cw.comment("rho = PW/(R*T)"),
    )

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckwt(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get molecular weight for all species"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT"
        + cc.sym
        + "( amrex::Real *  wt)",
    )
    cw.writer(fstream, "{")
    # call molecularWeight
    cw.writer(fstream, "get_mw(wt);")
    cw.writer(fstream, "}")


def ckmmwy(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("given y[species]: mass fractions"))
    cw.writer(fstream, cw.comment("s mean molecular weight (gm/mole)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWY"
        + cc.sym
        + "(amrex::Real *  y,  amrex::Real *  wtm)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "amrex::Real YOW = 0;")
    cw.writer(fstream, "amrex::Real tmp[%d];" % (nSpecies))
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)
    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "tmp[i] = y[i]*imw[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "YOW += tmp[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "")
    cw.writer(fstream, "*wtm = 1.0 / YOW;")
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckmmwx(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("given x[species]: mole fractions"))
    cw.writer(fstream, cw.comment("returns mean molecular weight (gm/mole)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWX"
        + cc.sym
        + "(amrex::Real *  x,  amrex::Real *  wtm)",
    )
    cw.writer(fstream, "{")
    cw.writer(
        fstream, "amrex::Real XW = 0;" + cw.comment(" see Eq 4 in CK Manual")
    )
    # molecular weights of all species
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "XW += x[%d]*%f; " % (spec_idx, species.weight)
            + cw.comment("%s" % species.name),
        )
    cw.writer(fstream, "*wtm = XW;")
    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckmmwc(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("given c[species]: molar concentration"))
    cw.writer(fstream, cw.comment("returns mean molecular weight (gm/mole)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWC"
        + cc.sym
        + "(amrex::Real *  c,  amrex::Real *  wtm)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(fstream, cw.comment("See Eq 5 in CK Manual"))
    cw.writer(fstream, "amrex::Real W = 0;")
    cw.writer(fstream, "amrex::Real sumC = 0;")
    # molecular weights of all species
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "W += c[%d]*%f; " % (spec_idx, species.weight)
            + cw.comment("%s" % species.name),
        )
    cw.writer(fstream)
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "sumC += c[id];")
    cw.writer(fstream, "}")
    cw.writer(
        fstream, cw.comment(" CK provides no guard against divison by zero")
    )
    cw.writer(fstream, "*wtm = W/sumC;")
    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckcpor(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get Cp/R as a function of T "))
    cw.writer(fstream, cw.comment("for all species (Eq 19)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPOR"
        + cc.sym
        + "(amrex::Real *  T, amrex::Real *  cpor)",
    )
    cw.writer(fstream, "{")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )

    # call routine
    cw.writer(fstream, "cp_R(cpor, tc);")
    cw.writer(fstream, "}")


def ckhort(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get H/RT as a function of T "))
    cw.writer(fstream, cw.comment("for all species (Eq 20)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHORT"
        + cc.sym
        + "(amrex::Real *  T, amrex::Real *  hort)",
    )
    cw.writer(fstream, "{")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )

    # call routine
    cw.writer(fstream, "speciesEnthalpy(hort, tc);")
    cw.writer(fstream, "}")


def cksor(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get S/R as a function of T "))
    cw.writer(fstream, cw.comment("for all species (Eq 21)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSOR"
        + cc.sym
        + "(amrex::Real *  T, amrex::Real *  sor)",
    )
    cw.writer(fstream, "{")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )

    # call routine
    cw.writer(fstream, "speciesEntropy(sor, tc);")
    cw.writer(fstream, "}")


def ckytx(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "convert y[species] (mass fracs) to x[species] (mole fracs)"
        ),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTX"
        + cc.sym
        + "(amrex::Real *  y,  amrex::Real *  x)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "amrex::Real YOW = 0;")
    cw.writer(fstream, "amrex::Real tmp[%d];" % (nSpecies))
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)
    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "tmp[i] = y[i]*imw[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "YOW += tmp[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "")
    cw.writer(fstream, "amrex::Real YOWINV = 1.0/YOW;")
    cw.writer(fstream, "")
    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "x[i] = y[i]*imw[i]*YOWINV;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckytcp(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "convert y[species] (mass fracs) to c[species] (molar conc)"
        ),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCP"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real"
        " *  c)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "amrex::Real YOW = 0;")
    cw.writer(fstream, "amrex::Real PWORT;")
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Compute inverse of mean molecular wt first")
    )
    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "c[i] = y[i]*imw[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "YOW += c[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "")
    cw.writer(fstream, cw.comment("PW/RT (see Eq. 7)"))
    cw.writer(
        fstream,
        "PWORT = (*P)/(YOW * %1.14e * (*T)); "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m,
    )

    # now compute conversion
    cw.writer(fstream, cw.comment("Now compute conversion"))
    cw.writer(fstream, "")
    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "c[i] = PWORT * y[i] * imw[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckytcr(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "convert y[species] (mass fracs) to c[species] (molar conc)"
        ),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCR"
        + cc.sym
        + "(amrex::Real *  rho, amrex::Real * /*T*/, amrex::Real * y, "
        " amrex::Real * c)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)
    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "c[i] = (*rho)  * y[i] * imw[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckxty(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "convert x[species] (mole fracs) to y[species] (mass fracs)"
        ),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTY"
        + cc.sym
        + "(amrex::Real *  x,  amrex::Real *  y)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "amrex::Real XW = 0; " + cw.comment("See Eq 4, 9 in CK Manual"),
    )

    # compute mean molecular weight first (eq 3)
    cw.writer(fstream, cw.comment("Compute mean molecular wt first"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "XW += x[%d]*%f; " % (spec_idx, species.weight)
            + cw.comment("%s" % species.name),
        )

    # now compute conversion
    cw.writer(fstream, cw.comment("Now compute conversion"))
    cw.writer(fstream, "amrex::Real XWinv = 1.0/XW;")
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "y[%d] = x[%d]*%f*XWinv; " % (spec_idx, spec_idx, species.weight),
        )

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckxtcp(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "convert x[species] (mole fracs) to c[species] (molar conc)"
        ),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCP"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real"
        " *  c)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(
        fstream,
        "amrex::Real PORT = (*P)/(%1.14e * (*T)); "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m
        + cw.comment("P/RT"),
    )
    # now compute conversion
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute conversion, see Eq 10"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "c[id] = x[id]*PORT;")
    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckxtcr(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "convert x[species] (mole fracs) to c[species] (molar conc)"
        ),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCR"
        + cc.sym
        + "(amrex::Real *  rho, amrex::Real * /*T*/, amrex::Real *  x,"
        " amrex::Real *  c)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(
        fstream,
        "amrex::Real XW = 0; " + cw.comment("See Eq 4, 11 in CK Manual"),
    )
    cw.writer(fstream, "amrex::Real ROW; ")

    # compute mean molecular weight first (eq 3)
    cw.writer(fstream, cw.comment("Compute mean molecular wt first"))
    for sp in range(nSpecies):
        species = species_info.nonqss_species[sp]
        cw.writer(
            fstream,
            "XW += x[%d]*%f; " % (species.idx, species.weight)
            + cw.comment("%s" % species.name),
        )

    # now compute conversion
    cw.writer(fstream, "ROW = (*rho) / XW;")
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute conversion, see Eq 11"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "c[id] = x[id]*ROW;")
    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckctx(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "convert c[species] (molar conc) to x[species] (mole fracs)"
        ),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTX"
        + cc.sym
        + "(amrex::Real *  c, amrex::Real *  x)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))
    cw.writer(fstream, "amrex::Real sumC = 0; ")

    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute sum of c "))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "sumC += c[id];")
    cw.writer(fstream, "}")

    # now compute conversion
    cw.writer(fstream)
    cw.writer(fstream, cw.comment(" See Eq 13 "))
    cw.writer(fstream, "amrex::Real sumCinv = 1.0/sumC;")
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "x[id] = c[id]*sumCinv;")
    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckcty(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "convert c[species] (molar conc) to y[species] (mass fracs)"
        ),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTY"
        + cc.sym
        + "(amrex::Real *  c, amrex::Real *  y)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream, "amrex::Real CW = 0; " + cw.comment("See Eq 12 in CK Manual")
    )

    # compute denominator in eq 12
    cw.writer(fstream, cw.comment("compute denominator in eq 12 first"))
    for sp in range(nSpecies):
        species = species_info.nonqss_species[sp]
        cw.writer(
            fstream,
            "CW += c[%d]*%f; " % (species.idx, species.weight)
            + cw.comment("%s" % species.name),
        )

    # now compute conversion
    cw.writer(fstream, cw.comment("Now compute conversion"))
    cw.writer(fstream, "amrex::Real CWinv = 1.0/CW;")
    for sp in range(nSpecies):
        species = species_info.nonqss_species[sp]
        cw.writer(
            fstream,
            "y[%d] = c[%d]*%f*CWinv; "
            % (species.idx, species.idx, species.weight),
        )

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ckcvml(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("get specific heat at constant volume as a function "),
    )
    cw.writer(fstream, cw.comment("of T for all species (molar units)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVML"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  cvml)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )

    # call routine
    cw.writer(fstream, "cv_R(cvml, tc);")

    # convert cv/R to cv
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(
        fstream,
        "cvml[id] *= %1.14e;"
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m,
    )
    cw.writer(fstream, "}")

    cw.writer(fstream, "}")


def ckcpml(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("get specific heat at constant pressure as a ")
    )
    cw.writer(
        fstream, cw.comment("function of T for all species (molar units)")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPML"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  cpml)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )

    # call routine
    cw.writer(fstream, "cp_R(cpml, tc);")

    # convert cp/R to cp
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(
        fstream,
        "cpml[id] *= %1.14e;"
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m,
    )
    cw.writer(fstream, "}")

    cw.writer(fstream, "}")


def ckuml(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get internal energy as a function "))
    cw.writer(fstream, cw.comment("of T for all species (molar units)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUML"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  uml)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m
        + cw.comment("R*T"),
    )

    # call routine
    cw.writer(fstream, "speciesInternalEnergy(uml, tc);")

    # convert e/RT to e with molar units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "uml[id] *= RT;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckhml(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("get enthalpy as a function "))
    cw.writer(fstream, cw.comment("of T for all species (molar units)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHML"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  hml)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m
        + cw.comment("R*T"),
    )

    # call routine
    cw.writer(fstream, "speciesEnthalpy(hml, tc);")

    # convert h/RT to h with molar units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "hml[id] *= RT;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckgml(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("get standard-state Gibbs energy as a function ")
    )
    cw.writer(fstream, cw.comment("of T for all species (molar units)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKGML"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  gml)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m
        + cw.comment("R*T"),
    )

    # call routine
    cw.writer(fstream, "gibbs(gml, tc);")

    # convert g/RT to g with molar units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "gml[id] *= RT;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckaml(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("get standard-state Helmholtz free energy as a ")
    )
    cw.writer(
        fstream, cw.comment("function of T for all species (molar units)")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKAML"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  aml)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m
        + cw.comment("R*T"),
    )

    # call routine
    cw.writer(fstream, "helmholtz(aml, tc);")

    # convert A/RT to A with molar units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "aml[id] *= RT;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def cksml(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("Returns the standard-state entropies in molar units"),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSML"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  sml)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )

    # call routine
    cw.writer(fstream, "speciesEntropy(sml, tc);")

    # convert s/R to s
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(
        fstream,
        "sml[id] *= %1.14e;"
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m,
    )
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckcvms(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the specific heats at constant volume")
    )
    cw.writer(fstream, cw.comment("in mass units (Eq. 29)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVMS"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  cvms)",
    )
    cw.writer(fstream, "{")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )

    # call routine
    cw.writer(fstream, "cv_R(cvms, tc);")

    # convert cv/R to cv with mass units
    cw.writer(fstream, cw.comment("multiply by R/molecularweight"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        ROW = (
            cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg
        ).m / species.weight
        cw.writer(
            fstream,
            "cvms[%d] *= %20.15e; " % (spec_idx, ROW)
            + cw.comment("%s" % species.name),
        )

    cw.writer(fstream, "}")


def ckcpms(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the specific heats at constant pressure")
    )
    cw.writer(fstream, cw.comment("in mass units (Eq. 26)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPMS"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  cpms)",
    )
    cw.writer(fstream, "{")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )

    # call routine
    cw.writer(fstream, "cp_R(cpms, tc);")

    # convert cp/R to cp with mass units
    cw.writer(fstream, cw.comment("multiply by R/molecularweight"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        ROW = (
            cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg
        ).m / species.weight
        cw.writer(
            fstream,
            "cpms[%d] *= %20.15e; " % (spec_idx, ROW)
            + cw.comment("%s" % species.name),
        )

    cw.writer(fstream, "}")


def ckums(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns internal energy in mass units (Eq 30.)")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUMS"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  ums)",
    )
    cw.writer(fstream, "{")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m
        + cw.comment("R*T"),
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # call routine
    cw.writer(fstream, "speciesInternalEnergy(ums, tc);")
    cw.writer(fstream)

    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "ums[i] *= RT*imw[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckhms(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns enthalpy in mass units (Eq 27.)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHMS"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  hms)",
    )
    cw.writer(fstream, "{")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; "
        + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m
        + cw.comment("R*T"),
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # call routine
    cw.writer(fstream, "speciesEnthalpy(hms, tc);")
    cw.writer(fstream)

    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "hms[i] *= RT*imw[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckgms(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns gibbs in mass units (Eq 31.)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKGMS"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  gms)",
    )
    cw.writer(fstream, "{")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m
        + cw.comment("R*T"),
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # call routine
    cw.writer(fstream, "gibbs(gms, tc);")
    cw.writer(fstream)

    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "gms[i] *= RT*imw[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckams(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns helmholtz in mass units (Eq 32.)"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKAMS"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  ams)",
    )
    cw.writer(fstream, "{")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )
    cw.writer(
        fstream,
        "amrex::Real RT = %1.14e*tT; "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m
        + cw.comment("R*T"),
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # call routine
    cw.writer(fstream, "helmholtz(ams, tc);")

    cw.writer(fstream, "for (int i = 0; i < %d; i++)" % (nSpecies))
    cw.writer(fstream, "{")
    cw.writer(fstream, "ams[i] *= RT*imw[i];")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def cksms(fstream, mechanism, species_info):
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the entropies in mass units (Eq 28.)")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSMS"
        + cc.sym
        + "(amrex::Real *  T,  amrex::Real *  sms)",
    )
    cw.writer(fstream, "{")

    # get temperature cache
    cw.writer(
        fstream, "amrex::Real tT = *T; " + cw.comment("temporary temperature")
    )
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT,"
        " tT*tT*tT*tT }; " + cw.comment("temperature cache"),
    )

    # call routine
    cw.writer(fstream, "speciesEntropy(sms, tc);")

    # convert s/R to s with mass units
    cw.writer(fstream, cw.comment("multiply by R/molecularweight"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        ROW = (
            cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg
        ).m / species.weight
        cw.writer(
            fstream,
            "sms[%d] *= %20.15e; " % (spec_idx, ROW)
            + cw.comment("%s" % species.name),
        )

    cw.writer(fstream, "}")


def ckwc(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("compute the production rate for each species")
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC"
        + cc.sym
        + "(amrex::Real *  T, amrex::Real *  C,  amrex::Real *  wdot)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))

    # convert C to SI units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to SI"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "C[id] *= 1.0e6;")
    cw.writer(fstream, "}")

    # call productionRate
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "productionRate(wdot, C, *T);")

    # convert C and wdot to chemkin units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "C[id] *= 1.0e-6;")
    cw.writer(fstream, "wdot[id] *= 1.0e-6;")
    cw.writer(fstream, "}")

    cw.writer(fstream, "}")


def ckwyp(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the molar production rate of species")
    )
    cw.writer(fstream, cw.comment("Given P, T, and mass fractions"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real"
        " *  wdot)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))

    cw.writer(
        fstream,
        "amrex::Real c[%d]; " % nSpecies + cw.comment("temporary storage"),
    )
    cw.writer(fstream, "amrex::Real YOW = 0; ")
    cw.writer(fstream, "amrex::Real PWORT; ")
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # compute inverse of mean molecular weight first (eq 3)
    cw.writer(
        fstream, cw.comment("Compute inverse of mean molecular wt first")
    )
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        species = species_info.nonqss_species[spec_idx]
        cw.writer(
            fstream,
            "YOW += y[%d]*imw[%d]; " % (spec_idx, spec_idx)
            + cw.comment("%s" % species.name),
        )

    cw.writer(fstream, cw.comment("PW/RT (see Eq. 7)"))
    cw.writer(
        fstream,
        "PWORT = (*P)/(YOW * %1.14e * (*T)); "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m,
    )

    cw.writer(fstream, cw.comment("multiply by 1e6 so c goes to SI"))
    cw.writer(fstream, "PWORT *= 1e6; ")

    # now compute conversion
    cw.writer(fstream, cw.comment("Now compute conversion (and go to SI)"))
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        cw.writer(
            fstream,
            "c[%d] = PWORT * y[%d]*imw[%d]; " % (spec_idx, spec_idx, spec_idx),
        )

    # call productionRate
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "productionRate(wdot, c, *T);")

    # convert wdot to chemkin units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "wdot[id] *= 1.0e-6;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckwxp(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the molar production rate of species")
    )
    cw.writer(fstream, cw.comment("Given P, T, and mole fractions"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP"
        + cc.sym
        + "(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real"
        " *  wdot)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))

    cw.writer(
        fstream,
        "amrex::Real c[%d]; " % nSpecies + cw.comment("temporary storage"),
    )

    cw.writer(
        fstream,
        "amrex::Real PORT = 1e6 * (*P)/(%1.14e * (*T)); "
        % ((cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg)).m
        + cw.comment("1e6 * P/RT so c goes to SI units"),
    )

    # now compute conversion
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute conversion, see Eq 10"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "c[id] = x[id]*PORT;")
    cw.writer(fstream, "}")

    # call productionRate
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "productionRate(wdot, c, *T);")

    # convert wdot to chemkin units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "wdot[id] *= 1.0e-6;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckwyr(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the molar production rate of species")
    )
    cw.writer(fstream, cw.comment("Given rho, T, and mass fractions"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR"
        + cc.sym
        + "(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, "
        " amrex::Real *  wdot)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))

    cw.writer(
        fstream,
        "amrex::Real c[%d]; " % nSpecies + cw.comment("temporary storage"),
    )
    cw.writer(fstream, "amrex::Real imw[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "get_imw(imw);")
    cw.writer(fstream)

    # now compute conversion
    cw.writer(
        fstream, cw.comment("See Eq 8 with an extra 1e6 so c goes to SI")
    )
    for sp in species_info.nonqss_species_list:
        spec_idx = species_info.ordered_idx_map[sp]
        cw.writer(
            fstream,
            "c[%d] = 1e6 * (*rho) * y[%d]*imw[%d]; "
            % (spec_idx, spec_idx, spec_idx),
        )

    # call productionRate
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("call productionRate"))
    cw.writer(fstream, "productionRate(wdot, c, *T);")

    # convert wdot to chemkin units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "wdot[id] *= 1.0e-6;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def ckwxr(fstream, mechanism, species_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the molar production rate of species")
    )
    cw.writer(fstream, cw.comment("Given rho, T, and mole fractions"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR"
        + cc.sym
        + "(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, "
        " amrex::Real *  wdot)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "int id; " + cw.comment("loop counter"))

    cw.writer(
        fstream,
        "amrex::Real c[%d]; " % nSpecies + cw.comment("temporary storage"),
    )

    cw.writer(
        fstream,
        "amrex::Real XW = 0; " + cw.comment("See Eq 4, 11 in CK Manual"),
    )
    cw.writer(fstream, "amrex::Real ROW; ")

    # compute mean molecular weight first (eq 3)
    cw.writer(fstream, cw.comment("Compute mean molecular wt first"))
    for sp in range(nSpecies):
        species = species_info.nonqss_species[sp]
        cw.writer(
            fstream,
            "XW += x[%d]*%f; " % (species.idx, species.weight)
            + cw.comment("%s" % species.name),
        )

    # now compute conversion
    cw.writer(fstream, cw.comment("Extra 1e6 factor to take c to SI"))
    cw.writer(fstream, "ROW = 1e6*(*rho) / XW;")
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Compute conversion, see Eq 11"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "c[id] = x[id]*ROW;")
    cw.writer(fstream, "}")

    # call productionRate
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "productionRate(wdot, c, *T);")

    # convert wdot to chemkin units
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("convert to chemkin units"))
    cw.writer(fstream, "for (id = 0; id < %d; ++id) {" % nSpecies)
    cw.writer(fstream, "wdot[id] *= 1.0e-6;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def T_given_ey(fstream):
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            " get temperature given internal energy in mass units and mass"
            " fracs"
        ),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void"
        " GET_T_GIVEN_EY(amrex::Real *  e, amrex::Real *  y, amrex::Real *  t,"
        " int * ierr)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "#ifdef CONVERGENCE")
    cw.writer(fstream, "const int maxiter = 5000;")
    cw.writer(fstream, "const amrex::Real tol  = 1.e-12;")
    cw.writer(fstream, "#else")
    cw.writer(fstream, "const int maxiter = 200;")
    cw.writer(fstream, "const amrex::Real tol  = 1.e-6;")
    cw.writer(fstream, "#endif")
    cw.writer(fstream, "amrex::Real ein  = *e;")
    cw.writer(
        fstream,
        "amrex::Real tmin = 90;"
        + cw.comment("max lower bound for thermo def"),
    )
    cw.writer(
        fstream,
        "amrex::Real tmax = 4000;"
        + cw.comment("min upper bound for thermo def"),
    )
    cw.writer(fstream, "amrex::Real e1,emin,emax,cv,t1,dt;")
    cw.writer(fstream, "int i;" + cw.comment(" loop counter"))
    cw.writer(fstream, "CKUBMS(&tmin, y, &emin);")
    cw.writer(fstream, "CKUBMS(&tmax, y, &emax);")
    cw.writer(fstream, "if (ein < emin) {")
    cw.writer(fstream, cw.comment("Linear Extrapolation below tmin"))
    cw.writer(fstream, "CKCVBS(&tmin, y, &cv);")
    cw.writer(fstream, "*t = tmin - (emin-ein)/cv;")
    cw.writer(fstream, "*ierr = 1;")
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "if (ein > emax) {")
    cw.writer(fstream, cw.comment("Linear Extrapolation above tmax"))
    cw.writer(fstream, "CKCVBS(&tmax, y, &cv);")
    cw.writer(fstream, "*t = tmax - (emax-ein)/cv;")
    cw.writer(fstream, "*ierr = 1;")
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "t1 = *t;")
    cw.writer(fstream, "if (t1 < tmin || t1 > tmax) {")
    cw.writer(fstream, "t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);")
    cw.writer(fstream, "}")
    cw.writer(fstream, "for (i = 0; i < maxiter; ++i) {")
    cw.writer(fstream, "CKUBMS(&t1,y,&e1);")
    cw.writer(fstream, "CKCVBS(&t1,y,&cv);")
    cw.writer(fstream, "dt = (ein - e1) / cv;")
    cw.writer(fstream, "if (dt > 100.) { dt = 100.; }")
    cw.writer(fstream, "else if (dt < -100.) { dt = -100.; }")
    cw.writer(fstream, "else if (fabs(dt) < tol) break;")
    cw.writer(fstream, "else if (t1+dt == t1) break;")
    cw.writer(fstream, "t1 += dt;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "*t = t1;")
    cw.writer(fstream, "*ierr = 0;")
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")
    cw.writer(fstream)


def T_given_hy(fstream):
    cw.writer(
        fstream,
        cw.comment(
            " get temperature given enthalpy in mass units and mass fracs"
        ),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void"
        " GET_T_GIVEN_HY(amrex::Real *  h, amrex::Real *  y, amrex::Real *  t,"
        " int * ierr)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "#ifdef CONVERGENCE")
    cw.writer(fstream, "const int maxiter = 5000;")
    cw.writer(fstream, "const amrex::Real tol  = 1.e-12;")
    cw.writer(fstream, "#else")
    cw.writer(fstream, "const int maxiter = 200;")
    cw.writer(fstream, "const amrex::Real tol  = 1.e-6;")
    cw.writer(fstream, "#endif")
    cw.writer(fstream, "amrex::Real hin  = *h;")
    cw.writer(
        fstream,
        "amrex::Real tmin = 90;"
        + cw.comment("max lower bound for thermo def"),
    )
    cw.writer(
        fstream,
        "amrex::Real tmax = 4000;"
        + cw.comment("min upper bound for thermo def"),
    )
    cw.writer(fstream, "amrex::Real h1,hmin,hmax,cp,t1,dt;")
    cw.writer(fstream, "int i;" + cw.comment(" loop counter"))
    cw.writer(fstream, "CKHBMS(&tmin, y, &hmin);")
    cw.writer(fstream, "CKHBMS(&tmax, y, &hmax);")
    cw.writer(fstream, "if (hin < hmin) {")
    cw.writer(fstream, cw.comment("Linear Extrapolation below tmin"))
    cw.writer(fstream, "CKCPBS(&tmin, y, &cp);")
    cw.writer(fstream, "*t = tmin - (hmin-hin)/cp;")
    cw.writer(fstream, "*ierr = 1;")
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "if (hin > hmax) {")
    cw.writer(fstream, cw.comment("Linear Extrapolation above tmax"))
    cw.writer(fstream, "CKCPBS(&tmax, y, &cp);")
    cw.writer(fstream, "*t = tmax - (hmax-hin)/cp;")
    cw.writer(fstream, "*ierr = 1;")
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "t1 = *t;")
    cw.writer(fstream, "if (t1 < tmin || t1 > tmax) {")
    cw.writer(fstream, "t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);")
    cw.writer(fstream, "}")
    cw.writer(fstream, "for (i = 0; i < maxiter; ++i) {")
    cw.writer(fstream, "CKHBMS(&t1,y,&h1);")
    cw.writer(fstream, "CKCPBS(&t1,y,&cp);")
    cw.writer(fstream, "dt = (hin - h1) / cp;")
    cw.writer(fstream, "if (dt > 100.) { dt = 100.; }")
    cw.writer(fstream, "else if (dt < -100.) { dt = -100.; }")
    cw.writer(fstream, "else if (fabs(dt) < tol) break;")
    cw.writer(fstream, "else if (t1+dt == t1) break;")
    cw.writer(fstream, "t1 += dt;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "*t = t1;")
    cw.writer(fstream, "*ierr = 0;")
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


# NEED TO DEAL WITH THIS WHEN QSS
def ckinu(fstream, mechanism, species_info, reaction_info):
    nReaction = mechanism.n_reactions

    maxsp = 0

    ns = [0 for _ in range(nReaction)]
    ki = [[] for _ in range(nReaction)]
    nu = [[] for _ in range(nReaction)]

    for orig_idx, _ in reaction_info.idxmap.items():
        reaction = mechanism.reaction(orig_idx)

        for symbol, coefficient in reaction.reactants.items():
            ki[orig_idx].append(species_info.ordered_idx_map[symbol])
            nu[orig_idx].append(-int(coefficient))
        for symbol, coefficient in reaction.products.items():
            ki[orig_idx].append(species_info.ordered_idx_map[symbol])
            nu[orig_idx].append(int(coefficient))

        maxsp = max(maxsp, len(ki[orig_idx]))

    for orig_idx, _ in reaction_info.idxmap.items():
        reaction = mechanism.reaction(orig_idx)

        ns[orig_idx] = len(ki[orig_idx])
        for _ in range(ns[orig_idx], maxsp):
            ki[orig_idx].append(0)
            nu[orig_idx].append(0)

    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "Returns a count of species in a reaction, and their indices"
        ),
    )
    cw.writer(fstream, cw.comment("and stoichiometric coefficients. (Eq 50)"))
    cw.writer(
        fstream,
        "void CKINU" + cc.sym + "(int * i, int * nspec, int * ki, int * nu)",
    )
    cw.writer(fstream, "{")

    str_ns = ",".join(str(x) for x in ns)
    cw.writer(fstream, "const int ns[%d] =\n     {%s};" % (nReaction, str_ns))

    str_ki = ",".join(
        ",".join(str(x) for x in ki[j]) for j in range(nReaction)
    )
    cw.writer(
        fstream,
        "const int kiv[%d] =\n     {%s};" % (nReaction * maxsp, str_ki),
    )

    str_nu = ",".join(
        ",".join(str(x) for x in nu[j]) for j in range(nReaction)
    )
    cw.writer(
        fstream,
        "const int nuv[%d] =\n     {%s};" % (nReaction * maxsp, str_nu),
    )

    cw.writer(fstream, "if (*i < 1) {")

    cw.writer(fstream, cw.comment("Return max num species per reaction"))
    cw.writer(fstream, "*nspec = %d;" % (maxsp))
    cw.writer(fstream, "} else {")
    cw.writer(fstream, "if (*i > %d) {" % (nReaction))
    cw.writer(fstream, "*nspec = -1;")
    cw.writer(fstream, "} else {")

    cw.writer(fstream, "*nspec = ns[*i-1];")
    cw.writer(fstream, "for (int j=0; j<*nspec; ++j) {")
    cw.writer(fstream, "ki[j] = kiv[(*i-1)*%d + j] + 1;" % maxsp)
    cw.writer(fstream, "nu[j] = nuv[(*i-1)*%d + j];" % maxsp)
    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")