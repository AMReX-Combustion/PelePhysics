"""Debug functions used for various subroutines."""
import re
import time

import symengine as sme

import ceptr.constants as cc
import ceptr.thermo as cth
import ceptr.writer as cw

# ##########
# Main debug functions
# ##########


def qssa_debug(
    fstream,
    mechanism,
    species_info,
    reaction_info,
    syms,
    helper_names_to_print,
    intermediate_names_to_print,
):
    """Run all of the QSSA debug routines."""

    print(f"Symbolic kf QSS print for debug")
    qssa_kf_debug(
        fstream,
        mechanism,
        species_info,
        reaction_info,
        syms,
    )

    print("Symbolic Sc qss print for debug")
    qssa_sc_qss_debug(
        fstream,
        mechanism,
        species_info,
        reaction_info,
        syms,
    )

    print("Symbolic qf qss print for debug")
    qssa_coeff_debug(
        fstream,
        mechanism,
        species_info,
        reaction_info,
        syms,
    )

    print("Symbolic qss terms print for debug")
    qssa_terms_debug(
        fstream,
        mechanism,
        species_info,
        reaction_info,
        syms,
        helper_names_to_print,
        intermediate_names_to_print,
    )

    print("Symbolic gibbs QSS print for debug")
    gibbsQSS_debug(
        fstream,
        mechanism,
        species_info,
        reaction_info,
        syms,
    )

    print("Symbolic enthalpy QSS print for debug")
    speciesEnthalpyQSS_debug(
        fstream,
        mechanism,
        species_info,
        reaction_info,
        syms,
    )


def thermo_debug(
    fstream,
    mechanism,
    species_info,
    reaction_info,
    syms,
):

    print("Symbolic gibbs print for debug")
    gibbs_debug(
        fstream,
        mechanism,
        species_info,
        reaction_info,
        syms,
    )

    print("Symbolic enthalpy print for debug")
    speciesEnthalpy_debug(
        fstream,
        mechanism,
        species_info,
        reaction_info,
        syms,
    )


def production_debug(fstream, mechanism, species_info, reaction_info, syms):

    print("Symbolic wdot print for debug")
    production_rate_debug(
        fstream,
        mechanism,
        species_info,
        reaction_info,
        syms,
    )


def jacobian_debug(
    fstream,
    mechanism,
    species_info,
    reaction_info,
    syms,
    dscqss_dscList=[],
    indexList=None,
):

    print("Symbolic dscqss_dsc term print for debug")
    dscqss_dsc_debug(
        fstream,
        mechanism,
        species_info,
        reaction_info,
        syms,
        dscqss_dscList,
        indexList,
    )

    print("Symbolic dscqss_dsc term print for debug")
    dscqss_dsc_fast_debug(
        fstream,
        mechanism,
        species_info,
        reaction_info,
        syms,
    )

    print("symbolic jacobian term print for debug")
    ajac_term_fast_debug(
        fstream,
        mechanism,
        species_info,
        reaction_info,
        syms,
    )


# ##########
# Sub-functions for debuging specific things
# ##########


def qssa_kf_debug(
    fstream,
    mechanism,
    species_info,
    reaction_info,
    syms,
):
    """Temporary QSSA kf debuging function."""
    n_species = species_info.n_species
    n_qssa_reactions = species_info.n_species
    cw.writer(fstream)
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_k_f_qss_debug"
        + "(const amrex::Real * tc, amrex::Real invT, amrex::Real * kf_qss)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "const amrex::Real T = tc[1];",
    )

    for ireac in range(reaction_info.n_qssa_reactions):

        cpp_str = syms.convert_to_cpp(syms.kf_qss_smp_tmp[ireac])
        cw.writer(
            fstream,
            "kf_qss[%s] = %s;"
            % (
                str(ireac),
                cpp_str,
            ),
        )

    cw.writer(fstream, "}")


def qssa_terms_debug(
    fstream,
    mechanism,
    species_info,
    reaction_info,
    syms,
    helper_names_to_print=[],
    intermediate_names_to_print=[],
):
    """Temporary QSSA terms debuging function."""
    n_species = species_info.n_species
    cw.writer(fstream)
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qss_terms_debug"
        + "(amrex::Real * sc, amrex::Real T)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream)
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T };"
        + cw.comment("temperature cache"),
    )
    cw.writer(fstream, "amrex::Real invT = 1.0 / tc[1];")
    cw.writer(fstream, "amrex::Real invT2 = invT * invT;")
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("reference concentration: P_atm / (RT) in inverse mol/m^3"),
    )
    cw.writer(
        fstream,
        "amrex::Real refC = %g / %g / T;"
        % (
            cc.Patm_pa,
            cc.R.to(cc.ureg.joule / (cc.ureg.mole / cc.ureg.kelvin)).m,
        ),
    )
    cw.writer(fstream, "amrex::Real refCinv = 1.0 / refC;")
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the mixture concentration"))
    cw.writer(fstream, "amrex::Real mixture = 0.0;")
    cw.writer(fstream, "for (int k = 0; k < %d; ++k) {" % n_species)
    cw.writer(fstream, "mixture += sc[k];")
    cw.writer(fstream, "}")
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
    cw.writer(fstream, "amrex::Real g_RT[%d];" % (n_species))
    cw.writer(fstream, "gibbs(g_RT, tc);")
    if species_info.n_qssa_species > 0:
        cw.writer(
            fstream,
            "amrex::Real g_RT_qss[%d];" % (species_info.n_qssa_species),
        )
        cw.writer(fstream, "gibbs_qss(g_RT_qss, tc);")
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the species enthalpy"))
    cw.writer(fstream, "amrex::Real h_RT[%d];" % (n_species))
    cw.writer(fstream, "speciesEnthalpy(h_RT, tc);")
    if species_info.n_qssa_species > 0:
        cw.writer(
            fstream,
            "amrex::Real h_RT_qss[%d];" % (species_info.n_qssa_species),
        )
        cw.writer(fstream, "speciesEnthalpy_qss(h_RT_qss, tc);")
    if species_info.n_qssa_species > 0:
        cw.writer(fstream)
        cw.writer(fstream, cw.comment("Fill sc_qss here"))
        # cw.writer(
        #    fstream, "amrex::Real sc_qss[%d];" % species_info.n_qssa_species
        # )
        cw.writer(
            fstream,
            "amrex::Real kf_qss[%d], qf_qss[%d], qr_qss[%d];"
            % (
                reaction_info.n_qssa_reactions,
                reaction_info.n_qssa_reactions,
                reaction_info.n_qssa_reactions,
            ),
        )
        cw.writer(fstream, "comp_k_f_qss(tc, invT, kf_qss);")
        cw.writer(
            fstream,
            "comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, tc, g_RT, g_RT_qss);",
        )

    for name in helper_names_to_print:
        cw.writer(
            fstream,
            "amrex::Real %s = %s;"
            % (
                name,
                sme.ccode(syms.intermediate_helpers_smp[name]),
            ),
        )
    for name in intermediate_names_to_print:
        cw.writer(
            fstream,
            "amrex::Real %s = %s;"
            % (
                name,
                sme.ccode(syms.intermediate_terms_smp[name]),
            ),
        )
    for name in helper_names_to_print:
        cw.writer(
            fstream,
            'std::cout << "  %s = " << %s << "\\n";' % (name, name),
        )
    for name in intermediate_names_to_print:
        cw.writer(
            fstream,
            'std::cout << "  %s = " << %s << "\\n";' % (name, name),
        )

    cw.writer(fstream, "}")


def qssa_coeff_debug(fstream, mechanism, species_info, reaction_info, syms):
    """Temporary QSSA coeff debuging function."""
    n_species = species_info.n_species
    cw.writer(fstream)
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qss_coeff_debug"
        + "(amrex::Real * sc, amrex::Real * qf_qss, amrex::Real * qr_qss, amrex::Real T)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream)
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T };"
        + cw.comment("temperature cache"),
    )
    cw.writer(fstream, "amrex::Real invT = 1.0 / tc[1];")
    cw.writer(fstream, "amrex::Real invT2 = invT * invT;")
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("reference concentration: P_atm / (RT) in inverse mol/m^3"),
    )
    cw.writer(
        fstream,
        "amrex::Real refC = %g / %g / T;"
        % (
            cc.Patm_pa,
            cc.R.to(cc.ureg.joule / (cc.ureg.mole / cc.ureg.kelvin)).m,
        ),
    )
    cw.writer(fstream, "amrex::Real refCinv = 1.0 / refC;")
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the mixture concentration"))
    cw.writer(fstream, "amrex::Real mixture = 0.0;")
    cw.writer(fstream, "for (int k = 0; k < %d; ++k) {" % n_species)
    cw.writer(fstream, "mixture += sc[k];")
    cw.writer(fstream, "}")
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
    cw.writer(fstream, "amrex::Real g_RT[%d];" % (n_species))
    cw.writer(fstream, "gibbs(g_RT, tc);")
    if species_info.n_qssa_species > 0:
        cw.writer(
            fstream,
            "amrex::Real g_RT_qss[%d];" % (species_info.n_qssa_species),
        )
        cw.writer(fstream, "gibbs_qss(g_RT_qss, tc);")
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the species enthalpy"))
    cw.writer(fstream, "amrex::Real h_RT[%d];" % (n_species))
    cw.writer(fstream, "speciesEnthalpy(h_RT, tc);")
    if species_info.n_qssa_species > 0:
        cw.writer(
            fstream,
            "amrex::Real h_RT_qss[%d];" % (species_info.n_qssa_species),
        )
        cw.writer(fstream, "speciesEnthalpy_qss(h_RT_qss, tc);")
    if species_info.n_qssa_species > 0:
        cw.writer(fstream)
        cw.writer(fstream, cw.comment("Fill sc_qss here"))
        # cw.writer(
        #    fstream, "amrex::Real sc_qss[%d];" % species_info.n_qssa_species
        # )
        cw.writer(
            fstream,
            "amrex::Real kf_qss[%d];" % (reaction_info.n_qssa_reactions,),
        )
        cw.writer(fstream, "comp_k_f_qss(tc, invT, kf_qss);")

    for ireac in range(reaction_info.n_qssa_reactions):
        times = time.time()
        cpp_str = syms.convert_to_cpp(syms.qf_qss_smp[ireac])
        timee = time.time()
        print("Made expr for qf %d (time = %.3g s)" % (ireac, timee - times))
        times = time.time()
        cw.writer(
            fstream,
            "qf_qss[%s] = %s;"
            % (
                str(ireac),
                cpp_str,
            ),
        )
        timee = time.time()
        print(
            "Printed expr for qf %d (time = %.3g s)" % (ireac, timee - times)
        )
    for ireac in range(reaction_info.n_qssa_reactions):
        times = time.time()
        cpp_str = syms.convert_to_cpp(syms.qr_qss_smp[ireac])
        timee = time.time()
        print("Made expr for qr %d (time = %.3g s)" % (ireac, timee - times))
        times = time.time()
        cw.writer(
            fstream,
            "qr_qss[%s] = %s;"
            % (
                str(ireac),
                cpp_str,
            ),
        )
        timee = time.time()
        print(
            "Printed expr for qr %d (time = %.3g s)" % (ireac, timee - times)
        )

    cw.writer(fstream, "}")


def qssa_sc_qss_debug(fstream, mechanism, species_info, reaction_info, syms):
    """Temporary SCQSS debuging function."""
    n_species = species_info.n_species
    cw.writer(fstream)
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_sc_qss_debug"
        + "(amrex::Real * sc, amrex::Real * sc_qss, amrex::Real T)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream)
    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T };"
        + cw.comment("temperature cache"),
    )
    cw.writer(fstream, "amrex::Real invT = 1.0 / tc[1];")
    cw.writer(fstream, "amrex::Real invT2 = invT * invT;")
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("reference concentration: P_atm / (RT) in inverse mol/m^3"),
    )
    cw.writer(
        fstream,
        "amrex::Real refC = %g / %g / T;"
        % (
            cc.Patm_pa,
            cc.R.to(cc.ureg.joule / (cc.ureg.mole / cc.ureg.kelvin)).m,
        ),
    )
    cw.writer(fstream, "amrex::Real refCinv = 1.0 / refC;")
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the mixture concentration"))
    cw.writer(fstream, "amrex::Real mixture = 0.0;")
    cw.writer(fstream, "for (int k = 0; k < %d; ++k) {" % n_species)
    cw.writer(fstream, "mixture += sc[k];")
    cw.writer(fstream, "}")
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
    cw.writer(fstream, "amrex::Real g_RT[%d];" % (n_species))
    cw.writer(fstream, "gibbs(g_RT, tc);")
    if species_info.n_qssa_species > 0:
        cw.writer(
            fstream,
            "amrex::Real g_RT_qss[%d];" % (species_info.n_qssa_species),
        )
        cw.writer(fstream, "gibbs_qss(g_RT_qss, tc);")
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the species enthalpy"))
    cw.writer(fstream, "amrex::Real h_RT[%d];" % (n_species))
    cw.writer(fstream, "speciesEnthalpy(h_RT, tc);")
    if species_info.n_qssa_species > 0:
        cw.writer(
            fstream,
            "amrex::Real h_RT_qss[%d];" % (species_info.n_qssa_species),
        )
        cw.writer(fstream, "speciesEnthalpy_qss(h_RT_qss, tc);")
    if species_info.n_qssa_species > 0:
        cw.writer(fstream)
        cw.writer(fstream, cw.comment("Fill sc_qss here"))
        # cw.writer(
        #    fstream, "amrex::Real sc_qss[%d];" % species_info.n_qssa_species
        # )
        cw.writer(
            fstream,
            "amrex::Real kf_qss[%d], qf_qss[%d], qr_qss[%d];"
            % (
                reaction_info.n_qssa_reactions,
                reaction_info.n_qssa_reactions,
                reaction_info.n_qssa_reactions,
            ),
        )
        cw.writer(fstream, "comp_k_f_qss(tc, invT, kf_qss);")
        cw.writer(
            fstream,
            "comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, tc, g_RT, g_RT_qss);",
        )

    # list_spec = [3,4,7,8,9,10,11,12,13,14,15,16,17]
    # list_spec = [10, 16, 17]
    # list_spec = list(range(species_info.n_qssa_species))
    list_spec = [0, 1, 2, 4, 5, 6, 7, 9, 11, 12, 13, 15, 8, 10, 14, 3, 16, 17]
    # list_spec = [1, 0]

    # for ispec in range(species_info.n_qssa_species):
    for ispec in list_spec:
        times = time.time()
        # Compute the common subexpressions using sympy
        sc_qss_cse = sme.cse([syms.sc_qss_smp[ispec]])

        # Write the reduced common expressions
        # The subexpressions are stored in cse index 0
        for cse_idx in range(len(sc_qss_cse[0])):
            common_exp = syms.convert_to_cpp(sc_qss_cse[0][cse_idx][1])
            common_exp = re.sub(
                r"(x)(\d{1,9})", r"x\2_" + str(ispec), common_exp
            )
            cw.writer(
                fstream,
                "const amrex::Real %s = %s;"
                % (
                    syms.convert_to_cpp(sc_qss_cse[0][cse_idx][0])
                    + "_"
                    + str(ispec),
                    common_exp,
                ),
            )

        # The full qss expression is stored in cse index 1
        cpp_str = syms.convert_to_cpp(sc_qss_cse[1])
        cpp_str = re.sub(r"(x)(\d{1,9})", r"x\2_" + str(ispec), cpp_str)
        timee = time.time()
        print("Made expr for spec %d (time = %.3g s)" % (ispec, timee - times))
        times = time.time()
        cw.writer(
            fstream,
            "sc_qss[%s] = %s;"
            % (
                str(ispec),
                cpp_str,
            ),
        )
        timee = time.time()
        print(
            "Printed expr for spec %d (time = %.3g s)" % (ispec, timee - times)
        )

    cw.writer(fstream, "}")


def gibbs_debug(fstream, mechanism, species_info, reaction_info, syms=None):
    """Temporary Write gibbs obtained with Sympy"""
    n_species = species_info.n_species
    n_reactions = mechanism.n_reactions

    cw.writer(fstream)

    # main function
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
        "gibbs_debug(amrex::Real * g_RT, const amrex::Real * tc)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "const amrex::Real T = tc[1];",
    )
    cw.writer(
        fstream,
        "const amrex::Real invT = 1.0 / tc[1];",
    )
    cw.writer(
        fstream,
        "const amrex::Real invT2 = invT * invT;",
    )

    for mid_temp in syms.midpointsList_sorted:
        cw.writer(
            fstream,
            "if (T < %.3g) {" % mid_temp,
        )
        for isymb, symb in enumerate(syms.g_RT_smp_tmp[mid_temp]["m"]):
            if str(symb).startswith("g_RT"):
                pass
            else:
                cpp_str = syms.convert_to_cpp(symb)
                cw.writer(
                    fstream,
                    "g_RT[%s] = %s;"
                    % (
                        str(isymb),
                        cpp_str,
                    ),
                )
        cw.writer(
            fstream,
            "}",
        )
        cw.writer(
            fstream,
            "else {",
        )
        for isymb, symb in enumerate(syms.g_RT_smp_tmp[mid_temp]["p"]):
            if str(symb).startswith("g_RT"):
                pass
            else:
                cpp_str = syms.convert_to_cpp(symb)
                cw.writer(
                    fstream,
                    "g_RT[%s] = %s;"
                    % (
                        str(isymb),
                        cpp_str,
                    ),
                )
        cw.writer(
            fstream,
            "}",
        )
    cw.writer(fstream, "}")


def gibbsQSS_debug(fstream, mechanism, species_info, reaction_info, syms=None):
    """Temporary Write gibbsQSS obtained with Sympy"""
    n_species = species_info.n_species
    n_reactions = mechanism.n_reactions

    cw.writer(fstream)

    # main function
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
        "gibbs_qss_debug(amrex::Real * g_RT_qss, const amrex::Real * tc)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "const amrex::Real T = tc[1];",
    )
    cw.writer(
        fstream,
        "const amrex::Real invT = 1.0 / tc[1];",
    )
    cw.writer(
        fstream,
        "const amrex::Real invT2 = invT * invT;",
    )

    for mid_temp in syms.midpointsQSSList_sorted:
        cw.writer(
            fstream,
            "if (T < %.3g) {" % mid_temp,
        )
        for isymb, symb in enumerate(syms.g_RT_qss_smp_tmp[mid_temp]["m"]):
            if str(symb).startswith("g_RT"):
                pass
            else:
                cpp_str = syms.convert_to_cpp(symb)
                cw.writer(
                    fstream,
                    "g_RT_qss[%s] = %s;"
                    % (
                        str(isymb),
                        cpp_str,
                    ),
                )
        cw.writer(
            fstream,
            "}",
        )
        cw.writer(
            fstream,
            "else {",
        )
        for isymb, symb in enumerate(syms.g_RT_qss_smp_tmp[mid_temp]["p"]):
            if str(symb).startswith("g_RT_qss"):
                pass
            else:
                cpp_str = syms.convert_to_cpp(symb)
                cw.writer(
                    fstream,
                    "g_RT_qss[%s] = %s;"
                    % (
                        str(isymb),
                        cpp_str,
                    ),
                )
        cw.writer(
            fstream,
            "}",
        )
    cw.writer(fstream, "}")


def speciesEnthalpy_debug(
    fstream, mechanism, species_info, reaction_info, syms=None
):
    """Temporary Write species enthalpy obtained with Sympy"""
    n_species = species_info.n_species
    n_reactions = mechanism.n_reactions

    cw.writer(fstream)

    # main function
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
        "speciesEnthalpy_debug(amrex::Real * h_RT, const amrex::Real * tc)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "const amrex::Real T = tc[1];",
    )
    cw.writer(
        fstream,
        "const amrex::Real invT = 1.0 / tc[1];",
    )
    cw.writer(
        fstream,
        "const amrex::Real invT2 = invT * invT;",
    )

    for mid_temp in syms.midpointsList_sorted:
        cw.writer(
            fstream,
            "if (T < %.3g) {" % mid_temp,
        )
        for isymb, symb in enumerate(syms.h_RT_smp_tmp[mid_temp]["m"]):
            if str(symb).startswith("h_RT"):
                pass
            else:
                cpp_str = syms.convert_to_cpp(symb)
                cw.writer(
                    fstream,
                    "h_RT[%s] = %s;"
                    % (
                        str(isymb),
                        cpp_str,
                    ),
                )
        cw.writer(
            fstream,
            "}",
        )
        cw.writer(
            fstream,
            "else {",
        )
        for isymb, symb in enumerate(syms.h_RT_smp_tmp[mid_temp]["p"]):
            if str(symb).startswith("h_RT"):
                pass
            else:
                cpp_str = syms.convert_to_cpp(symb)
                cw.writer(
                    fstream,
                    "h_RT[%s] = %s;"
                    % (
                        str(isymb),
                        cpp_str,
                    ),
                )
        cw.writer(
            fstream,
            "}",
        )
    cw.writer(fstream, "}")


def speciesEnthalpyQSS_debug(
    fstream, mechanism, species_info, reaction_info, syms=None
):
    """Temporary Write species enthalpy QSS obtained with Sympy"""
    n_species = species_info.n_species
    n_reactions = mechanism.n_reactions

    cw.writer(fstream)

    # main function
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
        "speciesEnthalpy_qss_debug(amrex::Real * h_RT_qss, const amrex::Real * tc)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "const amrex::Real T = tc[1];",
    )
    cw.writer(
        fstream,
        "const amrex::Real invT = 1.0 / tc[1];",
    )
    cw.writer(
        fstream,
        "const amrex::Real invT2 = invT * invT;",
    )

    for mid_temp in syms.midpointsQSSList_sorted:
        cw.writer(
            fstream,
            "if (T < %.3g) {" % mid_temp,
        )
        for isymb, symb in enumerate(syms.h_RT_qss_smp_tmp[mid_temp]["m"]):
            if str(symb).startswith("h_RT_qss"):
                pass
            else:
                cpp_str = syms.convert_to_cpp(symb)
                cw.writer(
                    fstream,
                    "h_RT_qss[%s] = %s;"
                    % (
                        str(isymb),
                        cpp_str,
                    ),
                )
        cw.writer(
            fstream,
            "}",
        )
        cw.writer(
            fstream,
            "else {",
        )
        for isymb, symb in enumerate(syms.h_RT_qss_smp_tmp[mid_temp]["p"]):
            if str(symb).startswith("h_RT_qss"):
                pass
            else:
                cpp_str = syms.convert_to_cpp(symb)
                cw.writer(
                    fstream,
                    "h_RT_qss[%s] = %s;"
                    % (
                        str(isymb),
                        cpp_str,
                    ),
                )
        cw.writer(
            fstream,
            "}",
        )
    cw.writer(fstream, "}")


def production_rate_debug(
    fstream, mechanism, species_info, reaction_info, syms=None
):
    """Temporary Write production rate obtained with Sympy. This is an expensive function"""
    n_species = species_info.n_species
    n_reactions = mechanism.n_reactions

    cw.writer(fstream)

    # main function
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
        " productionRate_debug(amrex::Real * wdot, amrex::Real * sc, amrex::Real T)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T };"
        + cw.comment("temperature cache"),
    )
    cw.writer(fstream, "const amrex::Real invT = 1.0 / tc[1];")
    cw.writer(fstream)

    if n_reactions == 0:
        cw.writer(fstream)
    else:
        cw.writer(
            fstream,
            cw.comment(
                "reference concentration: P_atm / (RT) in inverse mol/m^3"
            ),
        )
        cw.writer(
            fstream,
            "const amrex::Real refC = %g / %g * invT;"
            % (
                cc.Patm_pa,
                cc.R.to(cc.ureg.joule / (cc.ureg.mole / cc.ureg.kelvin)).m,
            ),
        )
        cw.writer(fstream, "const amrex::Real refCinv = 1 / refC;")

    if n_reactions > 0:
        # nclassd = n_reactions - nspecial
        # nCorr   = n3body + ntroe + nsri + nlindemann

        # Kc stuff
        cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
        cw.writer(fstream, "amrex::Real g_RT[%d];" % species_info.n_species)
        cw.writer(fstream, "gibbs(g_RT, tc);")
        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                "amrex::Real g_RT_qss[%d];" % (species_info.n_qssa_species),
            )
            cw.writer(fstream, "gibbs_qss(g_RT_qss, tc);")
        cw.writer(fstream)

        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                "amrex::Real sc_qss[%d];"
                % (max(1, species_info.n_qssa_species)),
            )
            cw.writer(
                fstream,
                "amrex::Real kf_qss[%d], qf_qss[%d], qr_qss[%d];"
                % (
                    reaction_info.n_qssa_reactions,
                    reaction_info.n_qssa_reactions,
                    reaction_info.n_qssa_reactions,
                ),
            )
            cw.writer(fstream, cw.comment("Fill sc_qss here"))
            cw.writer(fstream, "comp_k_f_qss(tc, invT, kf_qss);")
            # cw.writer(fstream,"comp_Kc_qss(invT, g_RT, g_RT_qss, Kc_qss);")
            cw.writer(
                fstream,
                "comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, tc, g_RT,"
                " g_RT_qss);",
            )
            cw.writer(fstream, "comp_sc_qss(sc_qss, qf_qss, qr_qss);")
            cw.writer(fstream)

            syms.write_array_to_cpp(syms.wdot_smp, "wdot", cw, fstream)

            cw.writer(fstream)

    cw.writer(fstream)

    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")

    cw.writer(fstream)


def ajac_term_debug(
    fstream,
    mechanism,
    species_info,
    reaction_info,
    syms=None,
    jacList=[],
    indexList=None,
):
    """Temporary Write jacobian term for debugging."""
    n_species = species_info.n_species
    n_reactions = mechanism.n_reactions

    cw.writer(fstream)

    # main
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
        " ajac_term_debug(amrex::Real * J, amrex::Real * sc, amrex::Real T, const int consP)",
    )
    cw.writer(fstream, "{")

    # Initialize the big Jacobian array
    cw.writer(fstream, "for (int i=0; i<%d; i++) {" % (n_species + 1) ** 2)
    cw.writer(fstream, "J[i] = 0.0;")
    cw.writer(fstream, "}")

    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T };"
        + cw.comment("temperature cache"),
    )
    cw.writer(fstream, "const amrex::Real invT = 1.0 / tc[1];")
    cw.writer(fstream)

    if n_reactions == 0:
        cw.writer(fstream)
    else:
        cw.writer(
            fstream,
            cw.comment(
                "reference concentration: P_atm / (RT) in inverse mol/m^3"
            ),
        )
        cw.writer(
            fstream,
            "const amrex::Real refC = %g / %g * invT;"
            % (
                cc.Patm_pa,
                cc.R.to(cc.ureg.joule / (cc.ureg.mole / cc.ureg.kelvin)).m,
            ),
        )
        cw.writer(fstream, "const amrex::Real refCinv = 1 / refC;")

    if n_reactions > 0:
        # nclassd = n_reactions - nspecial
        # nCorr   = n3body + ntroe + nsri + nlindemann

        # Kc stuff
        cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
        cw.writer(fstream, "amrex::Real g_RT[%d];" % species_info.n_species)
        cw.writer(fstream, "gibbs(g_RT, tc);")
        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                "amrex::Real g_RT_qss[%d];" % (species_info.n_qssa_species),
            )
            cw.writer(fstream, "gibbs_qss(g_RT_qss, tc);")
        cw.writer(fstream)
        cw.writer(fstream, cw.comment("compute the species enthalpy"))
        cw.writer(fstream, "amrex::Real h_RT[%d];" % (n_species))
        cw.writer(fstream, "speciesEnthalpy(h_RT, tc);")
        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                "amrex::Real h_RT_qss[%d];" % (species_info.n_qssa_species),
            )
            cw.writer(fstream, "speciesEnthalpy_qss(h_RT_qss, tc);")

        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                "amrex::Real sc_qss[%d];"
                % (max(1, species_info.n_qssa_species)),
            )
            cw.writer(
                fstream,
                "amrex::Real kf_qss[%d], qf_qss[%d], qr_qss[%d];"
                % (
                    reaction_info.n_qssa_reactions,
                    reaction_info.n_qssa_reactions,
                    reaction_info.n_qssa_reactions,
                ),
            )
            cw.writer(fstream, cw.comment("Fill sc_qss here"))
            cw.writer(fstream, "comp_k_f_qss(tc, invT, kf_qss);")
            # cw.writer(fstream,"comp_Kc_qss(invT, g_RT, g_RT_qss, Kc_qss);")
            cw.writer(
                fstream,
                "comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, tc, g_RT,"
                " g_RT_qss);",
            )
            cw.writer(fstream, "comp_sc_qss(sc_qss, qf_qss, qr_qss);")
            cw.writer(fstream)

            syms.write_array_to_cpp(jacList, f"J", cw, fstream, indexList)

            cw.writer(fstream)

    # dwdotdT
    cw.writer(fstream, "amrex::Real T_pert1, pertT;")
    cw.writer(
        fstream,
        "amrex::Real wdot_pert1[%d], wdot[%d];"
        % (
            n_species,
            n_species,
        ),
    )
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("dwdot/dT by finite difference"))
    cw.writer(fstream, "pertT = 1e-2;")
    cw.writer(fstream, "T_pert1 = T + pertT;")
    cw.writer(fstream)
    cw.writer(fstream, "productionRate(wdot_pert1, sc, T_pert1);")
    cw.writer(fstream, "productionRate(wdot, sc, T);")
    cw.writer(fstream)
    cw.writer(fstream, "for (int k = 0; k < %d ; k++) {" % n_species)
    cw.writer(
        fstream,
        "J[%d + k] = (wdot_pert1[k] - wdot[k])/(pertT);"
        % (n_species * (n_species + 1),),
    )
    cw.writer(fstream, "}")

    cw.writer(fstream)

    # depends on dwdotdT and dwdotdsc
    cw.writer(
        fstream,
        "amrex::Real c_R[%d], dcRdT[%d], e_RT[%d];"
        % (n_species, n_species, n_species),
    )
    cw.writer(fstream, "amrex::Real * eh_RT;")
    # if precond:
    #    cw.writer(fstream, "if (HP) {")
    # else:
    #    cw.writer(fstream, "if (consP) {")

    cw.writer(fstream, "if (consP) {")

    cw.writer(fstream, "cp_R(c_R, tc);")
    cw.writer(fstream, "dcvpRdT(dcRdT, tc);")
    cw.writer(fstream, "eh_RT = &h_RT[0];")

    cw.writer(fstream, "}")
    cw.writer(fstream, "else {")

    cw.writer(fstream, "cv_R(c_R, tc);")
    cw.writer(fstream, "dcvpRdT(dcRdT, tc);")
    cw.writer(fstream, "speciesInternalEnergy(e_RT, tc);")
    cw.writer(fstream, "eh_RT = &e_RT[0];")

    cw.writer(fstream, "}")

    cw.writer(fstream)

    cw.writer(
        fstream,
        "amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;",
    )
    cw.writer(fstream, "for (int k = 0; k < %d; ++k) {" % n_species)
    cw.writer(fstream, "cmix += c_R[k]*sc[k];")
    cw.writer(fstream, "dcmixdT += dcRdT[k]*sc[k];")
    cw.writer(fstream, "ehmix += eh_RT[k]*wdot[k];")
    cw.writer(
        fstream,
        "dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[%d+k];"
        % (n_species * (n_species + 1)),
    )
    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "amrex::Real cmixinv = 1.0/cmix;")
    cw.writer(fstream, "amrex::Real tmp1 = ehmix*cmixinv;")
    cw.writer(fstream, "amrex::Real tmp3 = cmixinv*T;")
    cw.writer(fstream, "amrex::Real tmp2 = tmp1*tmp3;")
    cw.writer(fstream, "amrex::Real dehmixdc;")

    cw.writer(fstream, cw.comment("dTdot/d[X]"))
    cw.writer(fstream, "for (int k = 0; k < %d; ++k) {" % n_species)
    cw.writer(fstream, "dehmixdc = 0.0;")
    cw.writer(fstream, "for (int m = 0; m < %d; ++m) {" % n_species)
    cw.writer(fstream, "dehmixdc += eh_RT[m]*J[k*%s+m];" % (n_species + 1))
    cw.writer(fstream, "}")
    cw.writer(
        fstream,
        "J[k*%d+%d] = tmp2*c_R[k] - tmp3*dehmixdc;"
        % (n_species + 1, n_species),
    )
    cw.writer(fstream, "}")

    cw.writer(fstream, cw.comment("dTdot/dT"))
    cw.writer(
        fstream,
        "J[%d] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;"
        % (n_species * (n_species + 1) + n_species),
    )

    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")

    cw.writer(fstream)


def dscqss_dsc_debug(
    fstream,
    mechanism,
    species_info,
    reaction_info,
    syms=None,
    dscqss_dscList=[],
    indexList=None,
):
    """Temporary Write dscqss_dsc for debugging."""
    n_species = species_info.n_species
    n_qssa_species = species_info.n_qssa_species
    n_reactions = mechanism.n_reactions

    cw.writer(fstream)

    # main
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
        " dscqss_dsc_debug(amrex::Real * dscqss_dsc, amrex::Real * sc, amrex::Real T)",
    )
    cw.writer(fstream, "{")

    # Initialize the big Jacobian array
    cw.writer(
        fstream, "for (int i=0; i<%d; i++) {" % (n_species * n_qssa_species)
    )
    cw.writer(fstream, "dscqss_dsc[i] = 0.0;")
    cw.writer(fstream, "}")

    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T };"
        + cw.comment("temperature cache"),
    )
    cw.writer(fstream, "const amrex::Real invT = 1.0 / tc[1];")
    cw.writer(fstream)

    if n_reactions == 0:
        cw.writer(fstream)
    else:
        cw.writer(
            fstream,
            cw.comment(
                "reference concentration: P_atm / (RT) in inverse mol/m^3"
            ),
        )
        cw.writer(
            fstream,
            "const amrex::Real refC = %g / %g * invT;"
            % (
                cc.Patm_pa,
                cc.R.to(cc.ureg.joule / (cc.ureg.mole / cc.ureg.kelvin)).m,
            ),
        )
        cw.writer(fstream, "const amrex::Real refCinv = 1 / refC;")

    if n_reactions > 0:
        # nclassd = n_reactions - nspecial
        # nCorr   = n3body + ntroe + nsri + nlindemann

        # Kc stuff
        cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
        cw.writer(fstream, "amrex::Real g_RT[%d];" % species_info.n_species)
        cw.writer(fstream, "gibbs(g_RT, tc);")
        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                "amrex::Real g_RT_qss[%d];" % (species_info.n_qssa_species),
            )
            cw.writer(fstream, "gibbs_qss(g_RT_qss, tc);")
        cw.writer(fstream)
        cw.writer(fstream, cw.comment("compute the species enthalpy"))
        cw.writer(fstream, "amrex::Real h_RT[%d];" % (n_species))
        cw.writer(fstream, "speciesEnthalpy(h_RT, tc);")
        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                "amrex::Real h_RT_qss[%d];" % (species_info.n_qssa_species),
            )
            cw.writer(fstream, "speciesEnthalpy_qss(h_RT_qss, tc);")

        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                "amrex::Real sc_qss[%d];"
                % (max(1, species_info.n_qssa_species)),
            )
            cw.writer(
                fstream,
                "amrex::Real kf_qss[%d], qf_qss[%d], qr_qss[%d];"
                % (
                    reaction_info.n_qssa_reactions,
                    reaction_info.n_qssa_reactions,
                    reaction_info.n_qssa_reactions,
                ),
            )
            cw.writer(fstream, cw.comment("Fill sc_qss here"))
            cw.writer(fstream, "comp_k_f_qss(tc, invT, kf_qss);")
            # cw.writer(fstream,"comp_Kc_qss(invT, g_RT, g_RT_qss, Kc_qss);")
            cw.writer(
                fstream,
                "comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, tc, g_RT,"
                " g_RT_qss);",
            )
            cw.writer(fstream, "comp_sc_qss(sc_qss, qf_qss, qr_qss);")
            cw.writer(fstream)

            syms.write_array_to_cpp(
                dscqss_dscList, f"dscqss_dsc", cw, fstream, indexList
            )

            cw.writer(fstream)

    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")

    cw.writer(fstream)


def dscqss_dsc_fast_debug(
    fstream,
    mechanism,
    species_info,
    reaction_info,
    syms=None,
):
    """Temporary Write dscqss_dsc for debugging."""

    n_species = species_info.n_species
    n_qssa_species = species_info.n_qssa_species
    n_reactions = mechanism.n_reactions

    cw.writer(fstream)

    # main
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
        " dscqss_dsc_fast_debug(amrex::Real * dscqss_dsc, amrex::Real * sc, amrex::Real T)",
    )
    cw.writer(fstream, "{")

    # Initialize the big Jacobian array
    cw.writer(
        fstream, "for (int i=0; i<%d; i++) {" % (n_species * n_qssa_species)
    )
    cw.writer(fstream, "dscqss_dsc[i] = 0.0;")
    cw.writer(fstream, "}")

    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T };"
        + cw.comment("temperature cache"),
    )
    cw.writer(fstream, "const amrex::Real invT = 1.0 / tc[1];")
    cw.writer(fstream)

    if n_reactions == 0:
        cw.writer(fstream)
    else:
        cw.writer(
            fstream,
            cw.comment(
                "reference concentration: P_atm / (RT) in inverse mol/m^3"
            ),
        )
        cw.writer(
            fstream,
            "const amrex::Real refC = %g / %g * invT;"
            % (
                cc.Patm_pa,
                cc.R.to(cc.ureg.joule / (cc.ureg.mole / cc.ureg.kelvin)).m,
            ),
        )
        cw.writer(fstream, "const amrex::Real refCinv = 1 / refC;")

    if n_reactions > 0:
        # nclassd = n_reactions - nspecial
        # nCorr   = n3body + ntroe + nsri + nlindemann

        # Kc stuff
        cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
        cw.writer(fstream, "amrex::Real g_RT[%d];" % species_info.n_species)
        cw.writer(fstream, "gibbs(g_RT, tc);")
        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                "amrex::Real g_RT_qss[%d];" % (species_info.n_qssa_species),
            )
            cw.writer(fstream, "gibbs_qss(g_RT_qss, tc);")
        cw.writer(fstream)
        cw.writer(fstream, cw.comment("compute the species enthalpy"))
        cw.writer(fstream, "amrex::Real h_RT[%d];" % (n_species))
        cw.writer(fstream, "speciesEnthalpy(h_RT, tc);")
        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                "amrex::Real h_RT_qss[%d];" % (species_info.n_qssa_species),
            )
            cw.writer(fstream, "speciesEnthalpy_qss(h_RT_qss, tc);")

        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                "amrex::Real sc_qss[%d];"
                % (max(1, species_info.n_qssa_species)),
            )
            cw.writer(
                fstream,
                "amrex::Real kf_qss[%d], qf_qss[%d], qr_qss[%d];"
                % (
                    reaction_info.n_qssa_reactions,
                    reaction_info.n_qssa_reactions,
                    reaction_info.n_qssa_reactions,
                ),
            )
            cw.writer(fstream, cw.comment("Fill sc_qss here"))
            cw.writer(fstream, "comp_k_f_qss(tc, invT, kf_qss);")
            # cw.writer(fstream,"comp_Kc_qss(invT, g_RT, g_RT_qss, Kc_qss);")
            cw.writer(
                fstream,
                "comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, tc, g_RT,"
                " g_RT_qss);",
            )
            cw.writer(fstream, "comp_sc_qss(sc_qss, qf_qss, qr_qss);")
            cw.writer(fstream)

            # Write the dscqss terms
            syms.write_dscqss_to_cpp(species_info, cw, fstream)

            cw.writer(fstream)

    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")

    cw.writer(fstream)
