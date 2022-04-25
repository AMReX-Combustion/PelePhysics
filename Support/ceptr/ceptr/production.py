import sys

import ceptr.constants as cc
import ceptr.utilities as cu
import ceptr.writer as cw


def productionRate(fstream, mechanism, species_info, reaction_info):
    nSpecies = species_info.nSpecies
    nReactions = mechanism.n_reactions

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

    if nReactions > 0:
        # nclassd = nReactions - nspecial
        # nCorr   = n3body + ntroe + nsri + nlindemann

        # reacs are sorted here
        for orig_idx, idx in reaction_info.idxmap.items():
            cw.writer(fstream)
            reaction = mechanism.reaction(orig_idx)
            cw.writer(
                fstream,
                cw.comment("reaction %d: %s" % (orig_idx, reaction.equation)),
            )
            # FIXME
            if hasattr(reaction, "ford"):
                # if len(reaction.ford) > 0:
                cw.writer(
                    fstream,
                    "qf[%d] = %s;"
                    % (
                        idx,
                        cu.QSSsortedPhaseSpace(
                            mechanism, species_info, reaction.ford
                        ),
                    ),
                )
            else:
                cw.writer(
                    fstream,
                    "qf[%d] = %s;"
                    % (
                        idx,
                        cu.QSSsortedPhaseSpace(
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
                        cu.QSSsortedPhaseSpace(
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
        cw.writer(fstream, "for (int i = 0; i < %d; ++i) {" % nSpecies)
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
        cw.writer(fstream, "amrex::Real g_RT[%d];" % species_info.nSpecies)
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
        cw.writer(fstream, "// Evaluate the kfs")
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

            KcExpArg = cu.sortedKcExpArg(mechanism, species_info, reaction)
            KcConvInv = cu.fKcConvInv(mechanism, species_info, reaction)

            dim = cu.phaseSpaceUnits(reaction.reactants)
            thirdBody = reaction.reaction_type == "three-body"
            falloff = reaction.reaction_type == "falloff"
            aeuc = cu.activationEnergyUnits()
            if not thirdBody and not falloff:
                uc = cu.prefactorUnits(
                    cc.ureg("mole/cm**3"), 1 - dim
                )  # Case 3 !PD, !TB
                ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), 1 - dim)
                A = (reaction.rate.pre_exponential_factor * ctuc).to(
                    uc.to_root_units()
                )
                beta = reaction.rate.temperature_exponent
                E = (
                    reaction.rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
            elif not falloff:
                uc = cu.prefactorUnits(
                    cc.ureg("mole/cm**3"), -dim
                )  # Case 2 !PD, TB
                ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), -dim)
                A = (reaction.rate.pre_exponential_factor * ctuc).to(
                    uc.to_root_units()
                )
                beta = reaction.rate.temperature_exponent
                E = (
                    reaction.rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
            else:
                uc = cu.prefactorUnits(
                    cc.ureg("mole/cm**3"), 1 - dim
                )  # Case 2 !PD, TB
                ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), 1 - dim)
                A = (reaction.high_rate.pre_exponential_factor * ctuc).to(
                    uc.to_root_units()
                )
                beta = reaction.high_rate.temperature_exponent
                E = (
                    reaction.high_rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)

                low_A = (reaction.low_rate.pre_exponential_factor * ctuc).to(
                    uc.to_root_units()
                )
                low_beta = reaction.low_rate.temperature_exponent
                low_E = (
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
                else:
                    print("Unrecognized reaction rate type")
                    sys.exit(1)

            cw.writer(
                fstream,
                "// (%d):  %s" % (orig_idx, reaction.equation),
            )
            cw.writer(fstream, "k_f = %.15g" % (A.m))
            if (beta == 0) and (E == 0):
                cw.writer(fstream, "           ;")
            else:
                if E == 0:
                    cw.writer(
                        fstream, "           * exp((%.15g) * tc[0]);" % (beta)
                    )
                elif beta == 0:
                    cw.writer(
                        fstream,
                        "           * exp(-(%.15g) * invT);"
                        % (((1.0 / cc.Rc / cc.ureg.kelvin)) * E),
                    )
                else:
                    cw.writer(
                        fstream,
                        "           * exp((%.15g) * tc[0] - (%.15g) * invT);"
                        % (beta, ((1.0 / cc.Rc / cc.ureg.kelvin)) * E),
                    )

            alpha = 1.0
            if not thirdBody and not falloff:
                # cw.writer(fstream,"Corr  = 1.0;")
                cw.writer(fstream, "qf[%d] *= k_f;" % idx)
            elif not falloff:
                alpha = enhancement_d_with_QSS(
                    mechanism, species_info, reaction
                )
                cw.writer(fstream, "Corr  = %s;" % (alpha))
                cw.writer(fstream, "qf[%d] *= Corr * k_f;" % idx)
            else:
                alpha = enhancement_d_with_QSS(
                    mechanism, species_info, reaction
                )
                cw.writer(fstream, "Corr  = %s;" % (alpha))
                cw.writer(
                    fstream,
                    "redP = Corr / k_f * 1e-%d * %.15g "
                    % (dim * 6, low_A.m * 10 ** (3**dim)),
                )
                cw.writer(
                    fstream,
                    "           * exp(%.15g  * tc[0] - %.15g  * (%.15g)"
                    " *invT);"
                    % (low_beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, low_E.m),
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
                elif reaction.sri:
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

            # FIXME
            if False:
                Ar, betar, Er = reaction.rev
                dim_rev = cu.phaseSpaceUnits(reaction.products)
                if not thirdBody and not falloff:
                    uc = cu.prefactorUnits(
                        cc.ureg("mole/cm**3"), 1 - dim_rev
                    )  # Case 3 !PD, !TB
                    ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), 1 - dim_rev)
                    print("Fixme grab rev params")
                elif not falloff:
                    uc = cu.prefactorUnits(
                        cc.ureg("mole/cm**3"), -dim_rev
                    )  # Case 2 !PD, TB
                    ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), -dim_rev)
                    print("Fixme grab rev params")
                else:
                    print("REV reaction cannot be PD")
                    sys.exit(1)
                cw.writer(fstream, "k_r = %.15g" % (Ar.m))
                if betar == 0:
                    cw.writer(
                        fstream,
                        "           * exp(- (%.15g) * invT);"
                        % ((1.0 / cc.Rc / cc.ureg.kelvin) * Er),
                    )
                else:
                    cw.writer(
                        fstream,
                        "           * exp(%.15g * tc[0] - (%.15g) * invT);"
                        % (betar, (1.0 / cc.Rc / cc.ureg.kelvin) * Er),
                    )
                if alpha == 1.0:
                    cw.writer(fstream, "qr[%d] *= k_r;" % idx)
                else:
                    cw.writer(fstream, "qr[%d] *= Corr * k_r;" % idx)
            else:
                if KcConvInv:
                    if alpha == 1.0:
                        cw.writer(
                            fstream,
                            "qr[%d] *= k_f * exp(-(%s)) * (%s);"
                            % (idx, KcExpArg, KcConvInv),
                        )
                    else:
                        cw.writer(
                            fstream,
                            "qr[%d] *= Corr * k_f * exp(-(%s)) * (%s);"
                            % (idx, KcExpArg, KcConvInv),
                        )
                else:
                    if alpha == 1.0:
                        cw.writer(
                            fstream,
                            "qr[%d] *= k_f * exp(-(%s));" % (idx, KcExpArg),
                        )
                    else:
                        cw.writer(
                            fstream,
                            "qr[%d] *= Corr * k_f * exp(-(%s));"
                            % (idx, KcExpArg),
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
        "const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T }; //"
        " temperature cache",
    )
    cw.writer(fstream, "const amrex::Real invT = 1.0 / tc[1];")
    cw.writer(fstream)

    if nReactions == 0:
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
    cw.writer(fstream, "for (int i = 0; i < %d; ++i) {" % nSpecies)
    cw.writer(fstream, "wdot[i] = 0.0;")
    cw.writer(fstream, "}")
    cw.writer(fstream)

    if nReactions > 0:
        # nclassd = nReactions - nspecial
        # nCorr   = n3body + ntroe + nsri + nlindemann

        # Mixt concentration for PD & TB
        cw.writer(fstream, cw.comment("compute the mixture concentration"))
        cw.writer(fstream, "amrex::Real mixture = 0.0;")
        cw.writer(fstream, "for (int i = 0; i < %d; ++i) {" % nSpecies)
        cw.writer(fstream, "mixture += sc[i];")
        cw.writer(fstream, "}")
        cw.writer(fstream)

        # Kc stuff
        cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
        cw.writer(fstream, "amrex::Real g_RT[%d];" % species_info.nSpecies)
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
            cw.writer(fstream, "// Fill sc_qss here")
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
            # FIXME
            if hasattr(reaction, "ford"):
                # if len(reaction.ford) > 0:
                forward_sc = cu.QSSsortedPhaseSpace(
                    mechanism, species_info, reaction.ford
                )
            else:
                forward_sc = cu.QSSsortedPhaseSpace(
                    mechanism, species_info, reaction.reactants
                )
            if reaction.reversible:
                reverse_sc = cu.QSSsortedPhaseSpace(
                    mechanism, species_info, reaction.products
                )
            else:
                reverse_sc = "0.0"

            KcExpArg = cu.sortedKcExpArg(mechanism, species_info, reaction)
            KcConvInv = cu.fKcConvInv(mechanism, species_info, reaction)

            dim = cu.phaseSpaceUnits(reaction.reactants)
            thirdBody = reaction.reaction_type == "three-body"
            falloff = reaction.reaction_type == "falloff"
            aeuc = cu.activationEnergyUnits()
            if not thirdBody and not falloff:
                uc = cu.prefactorUnits(
                    cc.ureg("mole/cm**3"), 1 - dim
                )  # Case 3 !PD, !TB
                ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), 1 - dim)
                A = (reaction.rate.pre_exponential_factor * ctuc).to(
                    uc.to_root_units()
                )
                beta = reaction.rate.temperature_exponent
                E = (
                    reaction.rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
            elif not falloff:
                uc = cu.prefactorUnits(
                    cc.ureg("mole/cm**3"), -dim
                )  # Case 2 !PD, TB
                ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), -dim)
                A = (reaction.rate.pre_exponential_factor * ctuc).to(
                    uc.to_root_units()
                )
                beta = reaction.rate.temperature_exponent
                E = (
                    reaction.rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
            else:
                uc = cu.prefactorUnits(
                    cc.ureg("mole/cm**3"), 1 - dim
                )  # Case 2 !PD, TB
                ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), 1 - dim)
                A = (reaction.high_rate.pre_exponential_factor * ctuc).to(
                    uc.to_root_units()
                )
                beta = reaction.high_rate.temperature_exponent
                E = (
                    reaction.high_rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)

                low_A = (reaction.low_rate.pre_exponential_factor * ctuc).to(
                    uc.to_root_units()
                )
                low_beta = reaction.low_rate.temperature_exponent
                low_E = (
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
                else:
                    print("Unrecognized reaction rate type")
                    sys.exit(1)

            cw.writer(
                fstream,
                "// (%d):  %s" % (orig_idx, reaction.equation),
            )
            cw.writer(fstream, "const amrex::Real k_f = %.15g" % (A.m))
            if (beta == 0) and (E == 0):
                cw.writer(fstream, "           ;")
            else:
                if E == 0:
                    cw.writer(
                        fstream, "           * exp((%.15g) * tc[0]);" % (beta)
                    )
                elif beta == 0:
                    cw.writer(
                        fstream,
                        "           * exp(-(%.15g) * invT);"
                        % (((1.0 / cc.Rc / cc.ureg.kelvin)) * E),
                    )
                else:
                    cw.writer(
                        fstream,
                        "           * exp((%.15g) * tc[0] - (%.15g) * invT);"
                        % (beta, (((1.0 / cc.Rc / cc.ureg.kelvin)) * E)),
                    )

            alpha = 1.0
            if not thirdBody and not falloff:
                cw.writer(
                    fstream,
                    "const amrex::Real qf = k_f * (%s);" % (forward_sc),
                )
            elif not falloff:
                alpha = enhancement_d_with_QSS(
                    mechanism, species_info, reaction
                )
                cw.writer(fstream, "const amrex::Real Corr = %s;" % (alpha))
                cw.writer(
                    fstream,
                    "const amrex::Real qf = Corr * k_f * (%s);" % (forward_sc),
                )
            else:
                alpha = enhancement_d_with_QSS(
                    mechanism, species_info, reaction
                )
                cw.writer(fstream, "amrex::Real Corr = %s;" % (alpha))
                cw.writer(
                    fstream,
                    "const amrex::Real redP = Corr / k_f * %.15g "
                    % (10 ** (-dim * 6) * low_A.m * 10 ** (3**dim)),
                )
                cw.writer(
                    fstream,
                    "           * exp(%.15g * tc[0] - %.15g * invT);"
                    % (low_beta, (1.0 / cc.Rc / cc.ureg.kelvin * low_E)),
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
                elif reaction.sri:
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

            # FIXME
            if False:
                Ar, betar, Er = reaction.rev
                dim_rev = cu.phaseSpaceUnits(reaction.products)
                if not thirdBody and not falloff:
                    uc = cu.prefactorUnits(
                        cc.ureg("mole/cm**3"), 1 - dim_rev
                    )  # Case 3 !PD, !TB
                    ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), 1 - dim_rev)
                    print("Fixme grab rev params")
                elif not falloff:
                    uc = cu.prefactorUnits(
                        cc.ureg("mole/cm**3"), -dim_rev
                    )  # Case 2 !PD, TB
                    ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), -dim_rev)
                    print("Fixme grab rev params")
                else:
                    print("REV reaction cannot be PD")
                    sys.exit(1)
                cw.writer(
                    fstream,
                    "const amrex::Real k_r = %.15g" % (Ar.m),
                )
                if betar == 0:
                    if Er == 0:
                        cw.writer(fstream, ";")
                    else:
                        cw.writer(
                            fstream,
                            "           * exp( - (%.15g) * invT);"
                            % ((1.0 / cc.Rc / cc.ureg.kelvin) * Er),
                        )
                else:
                    if Er == 0:
                        cw.writer(
                            fstream,
                            "           * exp(%.15g * tc[0]);" % (betar),
                        )
                    else:
                        cw.writer(
                            fstream,
                            "           * exp(%.15g * tc[0] - (%.15g) * invT);"
                            % (betar, (1.0 / cc.Rc / cc.ureg.kelvin) * Er),
                        )

                if alpha == 1.0:
                    cw.writer(
                        fstream,
                        "const amrex::Real qr = k_r * (%s);" % (reverse_sc),
                    )
                else:
                    cw.writer(
                        fstream,
                        "const amrex::Real qr = Corr * k_r * (%s);"
                        % (reverse_sc),
                    )
            else:
                if KcConvInv:
                    if alpha == 1.0:
                        cw.writer(
                            fstream,
                            "const amrex::Real qr = k_f * exp(-(%s)) * (%s) *"
                            " (%s);" % (KcExpArg, KcConvInv, reverse_sc),
                        )
                    else:
                        cw.writer(
                            fstream,
                            "const amrex::Real qr = Corr * k_f * exp(-(%s)) *"
                            " (%s) * (%s);"
                            % (KcExpArg, KcConvInv, reverse_sc),
                        )
                else:
                    if alpha == 1.0:
                        cw.writer(
                            fstream,
                            "const amrex::Real qr = k_f * exp(-(%s)) * (%s);"
                            % (KcExpArg, reverse_sc),
                        )
                    else:
                        cw.writer(
                            fstream,
                            "const amrex::Real qr = Corr * k_f * exp(-(%s)) *"
                            " (%s);" % (KcExpArg, reverse_sc),
                        )

            removeForward = cu.isRemoveForward(reaction_info, orig_idx)

            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(fstream, "//const amrex::Real qdot = qf - qr;")
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
            agents = sorted(
                agents,
                key=lambda x: next(
                    (y for y in species_info.all_species if y.name == x[0]),
                    None,
                ).idx,
            )
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


def enhancement_d_with_QSS(mechanism, species_info, reaction):
    thirdBody = reaction.reaction_type == "three-body"
    falloff = reaction.reaction_type == "falloff"
    if not thirdBody and not falloff:
        print("_enhancement_d called for a reaction without a third body")
        sys.exit(1)

    if not hasattr(reaction, "efficiencies"):
        print("FIXME")
        sys.exit(1)
        species, coefficient = thirdBody
        if species == "<mixture>":
            return "mixture"
        return "sc[%d]" % species_info.nonqss_species[species].idx

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
                species_info.ordered_idx_map[symbol] - species_info.nSpecies
            )
            if (efficiency - 1) == 1:
                alpha.append("%s" % (conc))
            else:
                alpha.append("%s*%s" % (factor, conc))

    return " + ".join(alpha).replace("+ -", "- ")
