"""Production functions."""
import sys

import ceptr.constants as cc
import ceptr.utilities as cu
import ceptr.writer as cw


def production_rate(fstream, mechanism, species_info, reaction_info):
    """Write production rate."""
    n_species = species_info.n_species
    n_reactions = mechanism.n_reactions

    if len(reaction_info.index) != 7:
        print("\n\nCheck this!!!\n")
        sys.exit(1)

    itroe = reaction_info.index[0:2]
    isri = reaction_info.index[1:3]
    ilindemann = reaction_info.index[2:4]
    # i3body = reaction_info.index[3:5]
    # isimple = reaction_info.index[4:6]
    # ispecial = reaction_info.index[5:7]

    ntroe = itroe[1] - itroe[0]
    nsri = isri[1] - isri[0]
    nlindemann = ilindemann[1] - ilindemann[0]
    # n3body = i3body[1] - i3body[0]
    # nsimple = isimple[1] - isimple[0]
    # nspecial = ispecial[1] - ispecial[0]

    # qdot
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "GPU version of productionRate: no more use of thermo namespace"
            " vectors "
        ),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qfqr(amrex::Real *"
        "  qf, amrex::Real * qr, amrex::Real * sc, amrex::Real * sc_qss,"
        " amrex::Real * tc, amrex::Real invT)",
    )
    cw.writer(fstream, "{")

    if n_reactions > 0:
        # nclassd = n_reactions - nspecial
        # nCorr   = n3body + ntroe + nsri + nlindemann

        # reacs are sorted here
        for orig_idx, idx in reaction_info.idxmap.items():
            cw.writer(fstream)
            reaction = mechanism.reaction(orig_idx)
            cw.writer(
                fstream,
                cw.comment("reaction %d: %s" % (orig_idx, reaction.equation)),
            )
            if bool(reaction.orders):
                cw.writer(
                    fstream,
                    "qf[%d] = %s;"
                    % (
                        idx,
                        cu.qss_sorted_phase_space(
                            mechanism, species_info, reaction.orders
                        ),
                    ),
                )
            else:
                cw.writer(
                    fstream,
                    "qf[%d] = %s;"
                    % (
                        idx,
                        cu.qss_sorted_phase_space(
                            mechanism, species_info, reaction.reactants
                        ),
                    ),
                )
            if reaction.reversible:
                cw.writer(
                    fstream,
                    "qr[%d] = %s;"
                    % (
                        idx,
                        cu.qss_sorted_phase_space(
                            mechanism, species_info, reaction.products
                        ),
                    ),
                )
            else:
                cw.writer(fstream, "qr[%d] = 0.0;" % (idx))

        cw.writer(fstream)

        # Mixt concentration for PD & TB
        cw.writer(fstream, cw.comment("compute the mixture concentration"))
        cw.writer(fstream, "amrex::Real mixture = 0.0;")
        cw.writer(fstream, "for (int i = 0; i < %d; ++i) {" % n_species)
        cw.writer(fstream, "mixture += sc[i];")
        cw.writer(fstream, "}")
        cw.writer(fstream)
        if species_info.nQSSspecies > 0:
            cw.writer(
                fstream,
                "for (int i = 0; i < %d; ++i) {" % species_info.nQSSspecies,
            )
            cw.writer(fstream, "mixture += sc_qss[i];")
            cw.writer(fstream, "}")
            cw.writer(fstream)

        # Kc stuff
        cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
        cw.writer(fstream, "amrex::Real g_RT[%d];" % species_info.n_species)
        cw.writer(fstream, "gibbs(g_RT, tc);")
        if species_info.nQSSspecies > 0:
            cw.writer(
                fstream,
                "amrex::Real g_RT_qss[%d];" % (species_info.nQSSspecies),
            )
            cw.writer(fstream, "gibbs_qss(g_RT_qss, tc);")

        cw.writer(fstream)

        cw.writer(
            fstream,
            cw.comment(
                "reference concentration: P_atm / (RT) in inverse mol/m^3"
            ),
        )
        cw.writer(
            fstream,
            "amrex::Real refC = %g / %g * invT;"
            % (
                cc.Patm_pa,
                cc.R.to(cc.ureg.joule / (cc.ureg.mole / cc.ureg.kelvin)).m,
            ),
        )
        cw.writer(fstream, "amrex::Real refCinv = 1 / refC;")

        cw.writer(fstream)

        # kfs
        cw.writer(fstream, cw.comment("Evaluate the kfs"))
        # cw.writer(fstream,"amrex::Real k_f[%d];"% nclassd)
        # cw.writer(fstream,"amrex::Real Corr[%d];" % nclassd)
        cw.writer(fstream, "amrex::Real k_f, k_r, Corr;")
        if ntroe > 0:
            cw.writer(
                fstream,
                "amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe,"
                " F_troe;",
            )
        if nsri > 0:
            cw.writer(fstream, "amrex::Real redP, F, X, F_sri;")
        cw.writer(fstream)

        # Loop like you're going through them in the mech.Linp order
        for orig_idx, idx in sorted(reaction_info.idxmap.items()):
            reaction = mechanism.reaction(orig_idx)

            kc_exp_arg = cu.sorted_kc_exp_arg(
                mechanism, species_info, reaction
            )
            kc_conv_inv = cu.fkc_conv_inv(mechanism, species_info, reaction)

            dim = cu.phase_space_units(reaction.reactants)
            third_body = reaction.reaction_type == "three-body"
            falloff = reaction.reaction_type == "falloff"
            is_troe = False
            is_sri = False
            aeuc = cu.activation_energy_units()
            if not third_body and not falloff:
                # Case 3 !PD, !TB
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
                pef = (
                    reaction.rate.pre_exponential_factor * ctuc
                ).to_base_units()
                beta = reaction.rate.temperature_exponent
                ae = (
                    reaction.rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
            elif not falloff:
                # Case 2 !PD, TB
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), -dim)
                pef = (
                    reaction.rate.pre_exponential_factor * ctuc
                ).to_base_units()
                beta = reaction.rate.temperature_exponent
                ae = (
                    reaction.rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
            else:
                # Case 2 !PD, TB
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
                pef = (
                    reaction.high_rate.pre_exponential_factor * ctuc
                ).to_base_units()
                beta = reaction.high_rate.temperature_exponent
                ae = (
                    reaction.high_rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)

                low_pef = (
                    reaction.low_rate.pre_exponential_factor * ctuc
                ).to_base_units()
                low_beta = reaction.low_rate.temperature_exponent
                low_ae = (
                    reaction.low_rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
                if reaction.rate.type == "Troe":
                    troe = reaction.rate.falloff_coeffs
                    ntroe = len(troe)
                    is_troe = True
                elif reaction.rate.type == "Sri":
                    sri = reaction.rate.falloff_coeffs
                    nsri = len(sri)
                    # is_sri = True
                elif reaction.rate.type == "Lindemann":
                    pass
                else:
                    print(
                        f"Unrecognized reaction rate type: {reaction.equation}"
                    )
                    sys.exit(1)

            cw.writer(
                fstream,
                cw.comment("(%d):  %s" % (orig_idx, reaction.equation)),
            )
            cw.writer(fstream, "k_f = %.15g" % (pef.m))
            if (beta == 0) and (ae == 0):
                cw.writer(fstream, "           ;")
            else:
                if ae == 0:
                    cw.writer(
                        fstream, "           * exp((%.15g) * tc[0]);" % (beta)
                    )
                elif beta == 0:
                    cw.writer(
                        fstream,
                        "           * exp(-(%.15g) * invT);"
                        % (((1.0 / cc.Rc / cc.ureg.kelvin)) * ae),
                    )
                else:
                    cw.writer(
                        fstream,
                        "           * exp((%.15g) * tc[0] - (%.15g) * invT);"
                        % (beta, ((1.0 / cc.Rc / cc.ureg.kelvin)) * ae),
                    )

            alpha = 1.0
            if not third_body and not falloff:
                # cw.writer(fstream,"Corr  = 1.0;")
                cw.writer(fstream, "qf[%d] *= k_f;" % idx)
            elif not falloff:
                alpha = enhancement_d_with_qss(
                    mechanism, species_info, reaction
                )
                cw.writer(fstream, "Corr  = %s;" % (alpha))
                cw.writer(fstream, "qf[%d] *= Corr * k_f;" % idx)
            else:
                alpha = enhancement_d_with_qss(
                    mechanism, species_info, reaction
                )
                cw.writer(fstream, "Corr  = %s;" % (alpha))
                cw.writer(
                    fstream,
                    "redP = Corr / k_f * 1e-%d * %.15g "
                    % (dim * 6, low_pef.m * 10 ** (3**dim)),
                )
                cw.writer(
                    fstream,
                    "           * exp(%.15g  * tc[0] - %.15g  * (%.15g)"
                    " *invT);"
                    % (low_beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, low_ae.m),
                )
                if is_troe:
                    cw.writer(fstream, "F = redP / (1.0 + redP);")
                    cw.writer(fstream, "logPred = log10(redP);")
                    cw.writer(fstream, "logFcent = log10(")
                    if abs(troe[1]) > 1.0e-100:
                        cw.writer(
                            fstream,
                            "    (%.15g)*exp(-tc[1] * %.15g)"
                            % (1.0 - troe[0], (1 / troe[1])),
                        )
                    else:
                        cw.writer(fstream, "     0.0 ")
                    if abs(troe[2]) > 1.0e-100:
                        cw.writer(
                            fstream,
                            "    + %.15g * exp(-tc[1] * %.15g)"
                            % (troe[0], (1 / troe[2])),
                        )
                    else:
                        cw.writer(fstream, "    + 0.0 ")
                    if ntroe == 4:
                        if troe[3] < 0:
                            cw.writer(
                                fstream, "    + exp(%.15g * invT));" % -troe[3]
                            )
                        else:
                            cw.writer(
                                fstream, "    + exp(-%.15g * invT));" % troe[3]
                            )
                    else:
                        cw.writer(fstream, "    + 0.0);")
                    cw.writer(fstream, "troe_c = -0.4 - 0.67 * logFcent;")
                    cw.writer(fstream, "troe_n = 0.75 - 1.27 * logFcent;")
                    cw.writer(
                        fstream,
                        "troe = (troe_c + logPred) / (troe_n - 0.14*(troe_c +"
                        " logPred));",
                    )
                    cw.writer(
                        fstream,
                        "F_troe = pow(10, logFcent / (1.0 + troe*troe));",
                    )
                    cw.writer(fstream, "Corr = F * F_troe;")
                    cw.writer(fstream, "qf[%d] *= Corr * k_f;" % idx)
                elif is_sri:
                    cw.writer(fstream, "F = redP / (1.0 + redP);")
                    cw.writer(fstream, "logPred = log10(redP);")
                    cw.writer(fstream, "X = 1.0 / (1.0 + logPred*logPred);")
                    if sri[1] < 0:
                        cw.writer(
                            fstream,
                            "F_sri = exp(X * log(%.15g * exp(%.15g*invT)"
                            % (sri[0], -sri[1]),
                        )
                    else:
                        cw.writer(
                            fstream,
                            "F_sri = exp(X * log(%.15g * exp(-%.15g*invT)"
                            % (sri[0], sri[1]),
                        )
                    if sri[2] > 1.0e-100:
                        cw.writer(fstream, "   +  exp(tc[0]/%.15g) " % sri[2])
                    else:
                        cw.writer(fstream, "   +  0. ")
                    cw.writer(
                        fstream,
                        "   *  (%d > 3 ? %.15g*exp(%.15g*tc[0]) : 1.0);"
                        % (nsri, sri[3], sri[4]),
                    )
                    cw.writer(fstream, "Corr = F * F_sri;")
                    cw.writer(fstream, "qf[%d] *= Corr * k_f;" % idx)
                elif nlindemann > 0:
                    cw.writer(fstream, "Corr = redP / (1. + redP);")
                    cw.writer(fstream, "qf[%d] *= Corr * k_f;" % idx)

            if kc_conv_inv:
                if alpha == 1.0:
                    cw.writer(
                        fstream,
                        "qr[%d] *= k_f * exp(-(%s)) * (%s);"
                        % (idx, kc_exp_arg, kc_conv_inv),
                    )
                else:
                    cw.writer(
                        fstream,
                        "qr[%d] *= Corr * k_f * exp(-(%s)) * (%s);"
                        % (idx, kc_exp_arg, kc_conv_inv),
                    )
            else:
                if alpha == 1.0:
                    cw.writer(
                        fstream,
                        "qr[%d] *= k_f * exp(-(%s));" % (idx, kc_exp_arg),
                    )
                else:
                    cw.writer(
                        fstream,
                        "qr[%d] *= Corr * k_f * exp(-(%s));"
                        % (idx, kc_exp_arg),
                    )

        cw.writer(fstream)

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")

    cw.writer(fstream)

    # main function
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
        " productionRate(amrex::Real * wdot, amrex::Real * sc, amrex::Real T)",
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

        if nsri > 0:
            cw.writer(fstream, "amrex::Real X, F_sri;")

    cw.writer(fstream)
    cw.writer(fstream, "for (int i = 0; i < %d; ++i) {" % n_species)
    cw.writer(fstream, "wdot[i] = 0.0;")
    cw.writer(fstream, "}")
    cw.writer(fstream)

    if n_reactions > 0:
        # nclassd = n_reactions - nspecial
        # nCorr   = n3body + ntroe + nsri + nlindemann

        # Mixt concentration for PD & TB
        cw.writer(fstream, cw.comment("compute the mixture concentration"))
        cw.writer(fstream, "amrex::Real mixture = 0.0;")
        cw.writer(fstream, "for (int i = 0; i < %d; ++i) {" % n_species)
        cw.writer(fstream, "mixture += sc[i];")
        cw.writer(fstream, "}")
        cw.writer(fstream)

        # Kc stuff
        cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
        cw.writer(fstream, "amrex::Real g_RT[%d];" % species_info.n_species)
        cw.writer(fstream, "gibbs(g_RT, tc);")
        if species_info.nQSSspecies > 0:
            cw.writer(
                fstream,
                "amrex::Real g_RT_qss[%d];" % (species_info.nQSSspecies),
            )
            cw.writer(fstream, "gibbs_qss(g_RT_qss, tc);")
        cw.writer(fstream)

        if species_info.nQSSspecies > 0:
            cw.writer(
                fstream,
                "amrex::Real sc_qss[%d];" % (max(1, species_info.nQSSspecies)),
            )
            cw.writer(
                fstream,
                "amrex::Real kf_qss[%d], qf_qss[%d], qr_qss[%d];"
                % (
                    species_info.nqssReactions,
                    species_info.nqssReactions,
                    species_info.nqssReactions,
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

        # Loop like you're going through them in the mech.Linp order
        for orig_idx, _ in reaction_info.idxmap.items():
            reaction = mechanism.reaction(orig_idx)
            cw.writer(fstream, "{")
            if bool(reaction.orders):
                forward_sc = cu.qss_sorted_phase_space(
                    mechanism, species_info, reaction.orders
                )
            else:
                forward_sc = cu.qss_sorted_phase_space(
                    mechanism, species_info, reaction.reactants
                )
            if reaction.reversible:
                reverse_sc = cu.qss_sorted_phase_space(
                    mechanism, species_info, reaction.products
                )
            else:
                reverse_sc = "0.0"

            kc_exp_arg = cu.sorted_kc_exp_arg(
                mechanism, species_info, reaction
            )
            kc_conv_inv = cu.fkc_conv_inv(mechanism, species_info, reaction)

            dim = cu.phase_space_units(reaction.reactants)
            third_body = reaction.reaction_type == "three-body"
            falloff = reaction.reaction_type == "falloff"
            is_troe = False
            is_sri = False
            aeuc = cu.activation_energy_units()
            if not third_body and not falloff:
                # Case 3 !PD, !TB
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
                pef = (
                    reaction.rate.pre_exponential_factor * ctuc
                ).to_base_units()
                beta = reaction.rate.temperature_exponent
                ae = (
                    reaction.rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
            elif not falloff:
                # Case 2 !PD, TB
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), -dim)
                pef = (
                    reaction.rate.pre_exponential_factor * ctuc
                ).to_base_units()
                beta = reaction.rate.temperature_exponent
                ae = (
                    reaction.rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
            else:
                # Case 2 !PD, TB
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
                pef = (
                    reaction.high_rate.pre_exponential_factor * ctuc
                ).to_base_units()
                beta = reaction.high_rate.temperature_exponent
                ae = (
                    reaction.high_rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)

                low_pef = (
                    reaction.low_rate.pre_exponential_factor * ctuc
                ).to_base_units()
                low_beta = reaction.low_rate.temperature_exponent
                low_ae = (
                    reaction.low_rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
                if reaction.rate.type == "Troe":
                    troe = reaction.rate.falloff_coeffs
                    ntroe = len(troe)
                    is_troe = True
                elif reaction.rate.type == "Sri":
                    sri = reaction.rate.falloff_coeffs
                    nsri = len(sri)
                    # is_sri = True
                elif reaction.rate.type == "Lindemann":
                    pass
                else:
                    print(
                        f"Unrecognized reaction rate type: {reaction.equation}"
                    )
                    sys.exit(1)

            cw.writer(
                fstream,
                cw.comment("(%d):  %s" % (orig_idx, reaction.equation)),
            )
            cw.writer(fstream, "const amrex::Real k_f = %.15g" % (pef.m))
            if (beta == 0) and (ae == 0):
                cw.writer(fstream, "           ;")
            else:
                if ae == 0:
                    cw.writer(
                        fstream, "           * exp((%.15g) * tc[0]);" % (beta)
                    )
                elif beta == 0:
                    cw.writer(
                        fstream,
                        "           * exp(-(%.15g) * invT);"
                        % (((1.0 / cc.Rc / cc.ureg.kelvin)) * ae),
                    )
                else:
                    cw.writer(
                        fstream,
                        "           * exp((%.15g) * tc[0] - (%.15g) * invT);"
                        % (beta, (((1.0 / cc.Rc / cc.ureg.kelvin)) * ae)),
                    )

            alpha = 1.0
            if not third_body and not falloff:
                cw.writer(
                    fstream,
                    "const amrex::Real qf = k_f * (%s);" % (forward_sc),
                )
            elif not falloff:
                alpha = enhancement_d_with_qss(
                    mechanism, species_info, reaction
                )
                cw.writer(fstream, "const amrex::Real Corr = %s;" % (alpha))
                cw.writer(
                    fstream,
                    "const amrex::Real qf = Corr * k_f * (%s);" % (forward_sc),
                )
            else:
                alpha = enhancement_d_with_qss(
                    mechanism, species_info, reaction
                )
                cw.writer(fstream, "amrex::Real Corr = %s;" % (alpha))
                cw.writer(
                    fstream,
                    "const amrex::Real redP = Corr / k_f * %.15g "
                    % (10 ** (-dim * 6) * low_pef.m * 10 ** (3**dim)),
                )
                cw.writer(
                    fstream,
                    "           * exp(%.15g * tc[0] - %.15g * invT);"
                    % (low_beta, (1.0 / cc.Rc / cc.ureg.kelvin * low_ae)),
                )
                if is_troe:
                    cw.writer(
                        fstream, "const amrex::Real F = redP / (1.0 + redP);"
                    )
                    cw.writer(
                        fstream, "const amrex::Real logPred = log10(redP);"
                    )
                    cw.writer(fstream, "const amrex::Real logFcent = log10(")
                    if abs(troe[1]) > 1.0e-100:
                        if 1.0 - troe[0] != 0:
                            cw.writer(
                                fstream,
                                "    %.15g * exp(-tc[1] * %.15g)"
                                % (1.0 - troe[0], (1 / troe[1])),
                            )
                    else:
                        cw.writer(fstream, "     0.0 ")
                    if abs(troe[2]) > 1.0e-100:
                        if troe[0] != 0:
                            cw.writer(
                                fstream,
                                "    + %.15g * exp(-tc[1] * %.15g)"
                                % (troe[0], (1 / troe[2])),
                            )
                    else:
                        cw.writer(fstream, "    + 0.0 ")
                    if ntroe == 4:
                        if troe[3] < 0:
                            cw.writer(
                                fstream, "    + exp(%.15g * invT));" % -troe[3]
                            )
                        else:
                            cw.writer(
                                fstream, "    + exp(-%.15g * invT));" % troe[3]
                            )
                    else:
                        cw.writer(fstream, "    + 0.0);")
                    cw.writer(
                        fstream,
                        "const amrex::Real troe_c = -0.4 - 0.67 * logFcent;",
                    )
                    cw.writer(
                        fstream,
                        "const amrex::Real troe_n = 0.75 - 1.27 * logFcent;",
                    )
                    cw.writer(
                        fstream,
                        "const amrex::Real troe = (troe_c + logPred) / (troe_n"
                        " - 0.14 * (troe_c + logPred));",
                    )
                    cw.writer(
                        fstream,
                        "const amrex::Real F_troe = pow(10, logFcent / (1.0 +"
                        " troe * troe));",
                    )
                    cw.writer(fstream, "Corr = F * F_troe;")
                    cw.writer(
                        fstream,
                        "const amrex::Real qf = Corr * k_f * (%s);"
                        % (forward_sc),
                    )
                elif is_sri:
                    cw.writer(
                        fstream, "const amrex::Real F = redP / (1.0 + redP);"
                    )
                    cw.writer(
                        fstream, "const amrex::Real logPred = log10(redP);"
                    )
                    cw.writer(fstream, "X = 1.0 / (1.0 + logPred*logPred);")
                    if sri[1] < 0:
                        cw.writer(
                            fstream,
                            "F_sri = exp(X * log(%.15g * exp(%.15g * invT)"
                            % (sri[0], -sri[1]),
                        )
                    else:
                        cw.writer(
                            fstream,
                            "F_sri = exp(X * log(%.15g * exp(-%.15g * invT)"
                            % (sri[0], sri[1]),
                        )
                    if sri[2] > 1.0e-100:
                        cw.writer(
                            fstream, "   +  exp(tc[0] / %.15g) " % sri[2]
                        )
                    else:
                        cw.writer(fstream, "   +  0. ")
                    cw.writer(
                        fstream,
                        "   *  (%d > 3 ? %.15g * exp(%.15g * tc[0]) : 1.0);"
                        % (nsri, sri[3], sri[4]),
                    )
                    cw.writer(fstream, "Corr = F * F_sri;")
                    cw.writer(
                        fstream,
                        "const amrex::Real qf = Corr * k_f * (%s);"
                        % (forward_sc),
                    )
                elif nlindemann > 0:
                    cw.writer(fstream, "Corr = redP / (1.0 + redP);")
                    cw.writer(
                        fstream,
                        "const amrex::Real qf = Corr * k_f * (%s);"
                        % (forward_sc),
                    )

            if kc_conv_inv:
                if alpha == 1.0:
                    cw.writer(
                        fstream,
                        "const amrex::Real qr = k_f * exp(-(%s)) * (%s) *"
                        " (%s);" % (kc_exp_arg, kc_conv_inv, reverse_sc),
                    )
                else:
                    cw.writer(
                        fstream,
                        "const amrex::Real qr = Corr * k_f * exp(-(%s)) *"
                        " (%s) * (%s);"
                        % (kc_exp_arg, kc_conv_inv, reverse_sc),
                    )
            else:
                if alpha == 1.0:
                    cw.writer(
                        fstream,
                        "const amrex::Real qr = k_f * exp(-(%s)) * (%s);"
                        % (kc_exp_arg, reverse_sc),
                    )
                else:
                    cw.writer(
                        fstream,
                        "const amrex::Real qr = Corr * k_f * exp(-(%s)) *"
                        " (%s);" % (kc_exp_arg, reverse_sc),
                    )

            remove_forward = cu.is_remove_forward(reaction_info, orig_idx)

            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(
                    fstream, cw.comment("const amrex::Real qdot = qf - qr;")
                )
                cw.writer(fstream, "const amrex::Real qdot = - qr;")
            else:
                cw.writer(fstream, "const amrex::Real qdot = qf - qr;")

            reaction = mechanism.reaction(orig_idx)
            lst_reactants = [(k, v) for k, v in reaction.reactants.items()]
            lst_products = [(k, v) for k, v in reaction.products.items()]
            all_agents = list(set(lst_reactants + lst_products))
            agents = []
            # remove QSS species from agents
            for symbol, coefficient in all_agents:
                if symbol not in species_info.qss_species_list:
                    agents.append((symbol, coefficient))
            dict_species = {
                v: i for i, v in enumerate(species_info.all_species_list)
            }
            agents = sorted(agents, key=lambda v: dict_species[v[0]])
            # note that a species might appear as both reactant and product
            # a species might alos appear twice or more on on each side
            # agents is a set that contains unique (symbol, coefficient)
            for a in agents:
                symbol, coefficient = a
                for b in reaction.reactants:
                    if b == a[0]:
                        if coefficient == 1.0:
                            cw.writer(
                                fstream,
                                "wdot[%d] -= qdot;"
                                % (species_info.ordered_idx_map[symbol]),
                            )
                        else:
                            cw.writer(
                                fstream,
                                "wdot[%d] -= %f * qdot;"
                                % (
                                    species_info.ordered_idx_map[symbol],
                                    coefficient,
                                ),
                            )
                for b in reaction.products:
                    if b == a[0]:
                        if coefficient == 1.0:
                            cw.writer(
                                fstream,
                                "wdot[%d] += qdot;"
                                % (species_info.ordered_idx_map[symbol]),
                            )
                        else:
                            cw.writer(
                                fstream,
                                "wdot[%d] += %f * qdot;"
                                % (
                                    species_info.ordered_idx_map[symbol],
                                    coefficient,
                                ),
                            )
            cw.writer(fstream, "}")
            cw.writer(fstream)

    cw.writer(fstream)

    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")

    cw.writer(fstream)


def enhancement_d_with_qss(mechanism, species_info, reaction):
    """Write enhancement with QSS."""
    third_body = reaction.reaction_type == "three-body"
    falloff = reaction.reaction_type == "falloff"
    if not third_body and not falloff:
        print("enhancement_d called for a reaction without a third body")
        sys.exit(1)

    if not hasattr(reaction, "efficiencies"):
        print("FIXME EFFICIENCIES")
        sys.exit(1)
        species, coefficient = third_body
        if species == "<mixture>":
            return "mixture"
        return "sc[%d]" % species_info.ordered_idx_map[species]

    efficiencies = reaction.efficiencies
    alpha = ["mixture"]
    for _, (symbol, efficiency) in enumerate(efficiencies.items()):
        if symbol not in species_info.qss_species_list:
            factor = "(%.15g)" % (efficiency - 1)
            conc = "sc[%d]" % species_info.ordered_idx_map[symbol]
            if (efficiency - 1) == 1:
                alpha.append("%s" % (conc))
            else:
                alpha.append("%s*%s" % (factor, conc))
        else:
            factor = "(%.15g)" % (efficiency - 1)
            conc = "sc_qss[%d]" % (
                species_info.ordered_idx_map[symbol] - species_info.n_species
            )
            if (efficiency - 1) == 1:
                alpha.append("%s" % (conc))
            else:
                alpha.append("%s*%s" % (factor, conc))

    return " + ".join(alpha).replace("+ -", "- ")
