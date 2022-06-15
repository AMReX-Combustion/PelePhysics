"""Write jacobian functions."""
import copy
import sys
from collections import Counter, OrderedDict

import ceptr.constants as cc
import ceptr.utilities as cu
import ceptr.writer as cw


def ajac(
    fstream, mechanism, species_info, reaction_info, precond=False, syms=None
):
    """Write jacobian for a reaction."""
    n_species = species_info.n_species

    cw.writer(fstream)
    if precond:
        cw.writer(
            fstream, cw.comment("compute an approx to the reaction Jacobian")
        )
    else:
        cw.writer(fstream, cw.comment("compute the reaction Jacobian"))
    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
    if precond:
        cw.writer(
            fstream,
            "void aJacobian_precond(amrex::Real *  J, amrex::Real *  sc, amrex::Real"
            " T, const int HP)",
        )
    else:
        cw.writer(
            fstream,
            "void aJacobian(amrex::Real * J, amrex::Real * sc, amrex::Real T,"
            " const int consP)",
        )
    cw.writer(fstream, "{")

    cw.writer(fstream)
    # Analytical jacobian not ready with QSS
    if species_info.n_qssa_species > 0:
        cw.writer(
            fstream, cw.comment("Do not use Analytical Jacobian with QSSA")
        )
        cw.writer(fstream, "//amrex::Abort();")
        cw.writer(fstream)

    cw.writer(fstream, "for (int i=0; i<%d; i++) {" % (n_species + 1) ** 2)
    cw.writer(fstream, "J[i] = 0.0;")
    cw.writer(fstream, "}")

    cw.writer(fstream)

    cw.writer(fstream, "amrex::Real wdot[%d];" % (n_species))
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (n_species))
    cw.writer(fstream, "wdot[k] = 0.0;")
    cw.writer(fstream, "}")

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
        cw.writer(
            fstream, "amrex::Real sc_qss[%d];" % species_info.n_qssa_species
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
        cw.writer(fstream, "comp_k_f_qss(tc, invT, kf_qss);")
        cw.writer(
            fstream,
            "comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, tc, g_RT, g_RT_qss);",
        )
        cw.writer(fstream, "comp_sc_qss(sc_qss, qf_qss, qr_qss);")
        cw.writer(fstream)

    cw.writer(fstream)

    cw.writer(
        fstream,
        "amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;",
    )
    cw.writer(fstream, "amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;")
    cw.writer(fstream, "amrex::Real dqdci, dcdc_fac, dqdc[%d];" % (n_species))
    cw.writer(fstream, "amrex::Real Pr, fPr, F, k_0, logPr;")
    cw.writer(
        fstream,
        "amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;",
    )
    cw.writer(fstream, "amrex::Real Fcent1, Fcent2, Fcent3, Fcent;")
    cw.writer(fstream, "amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;")
    cw.writer(
        fstream,
        "amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr,"
        " dlnCorrdT;",
    )
    cw.writer(fstream, "const amrex::Real ln10 = log(10.0);")
    cw.writer(fstream, "const amrex::Real log10e = 1.0/log(10.0);")

    for orig_idx, _ in reaction_info.idxmap.items():
        # if orig_idx == 35:
        #     exit()
        reaction = mechanism.reaction(orig_idx)

        cw.writer(
            fstream,
            cw.comment("reaction %d: %s" % (orig_idx, reaction.equation)),
        )
        ajac_reaction_d(
            fstream,
            mechanism,
            species_info,
            reaction_info,
            reaction,
            orig_idx,
            precond=precond,
            syms=syms,
        )
        cw.writer(fstream)

    cw.writer(
        fstream,
        "amrex::Real c_R[%d], dcRdT[%d], e_RT[%d];"
        % (n_species, n_species, n_species),
    )
    cw.writer(fstream, "amrex::Real * eh_RT;")
    if precond:
        cw.writer(fstream, "if (HP) {")
    else:
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

    cw.writer(fstream, "}")


def ajac_reaction_d(
    fstream,
    mechanism,
    species_info,
    reaction_info,
    reaction,
    orig_idx,
    precond=False,
    syms=None,
):
    """Write jacobian of reaction."""
    n_species = species_info.n_species
    remove_forward = cu.is_remove_forward(reaction_info, orig_idx)

    dim = cu.phase_space_units(reaction.reactants)
    third_body = reaction.reaction_type == "three-body"
    falloff = reaction.reaction_type == "falloff"
    is_sri = False
    is_troe = False
    aeuc = cu.activation_energy_units()
    if not third_body and not falloff:
        # Case 3 !PD, !TB
        cw.writer(
            fstream,
            cw.comment("a non-third-body and non-pressure-fall-off reaction"),
        )
        ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
        pef = (reaction.rate.pre_exponential_factor * ctuc).to_base_units()
        beta = reaction.rate.temperature_exponent
        ae = (
            reaction.rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
        ).to(aeuc)
    elif not falloff:
        # Case 2 !PD, TB
        cw.writer(
            fstream,
            cw.comment("a third-body and non-pressure-fall-off reaction"),
        )
        ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), -dim)
        pef = (reaction.rate.pre_exponential_factor * ctuc).to_base_units()
        beta = reaction.rate.temperature_exponent
        ae = (
            reaction.rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
        ).to(aeuc)
    else:
        # Case 1 PD, TB
        cw.writer(
            fstream,
            cw.comment("a third-body and pressure-fall-off reaction"),
        )
        ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
        pef = (
            reaction.rate.high_rate.pre_exponential_factor * ctuc
        ).to_base_units()
        beta = reaction.rate.high_rate.temperature_exponent
        ae = (
            reaction.rate.high_rate.activation_energy
            * cc.ureg.joule
            / cc.ureg.kmol
        ).to(aeuc)

        low_pef = (
            reaction.rate.low_rate.pre_exponential_factor * ctuc
        ).to_base_units()
        low_beta = reaction.rate.low_rate.temperature_exponent
        low_ae = (
            reaction.rate.low_rate.activation_energy
            * cc.ureg.joule
            / cc.ureg.kmol
        ).to(aeuc)
        if reaction.rate.type == "Troe":
            troe = reaction.rate.falloff_coeffs
            ntroe = len(troe)
            is_troe = True
        elif reaction.rate.type == "Sri":
            pass
            # sri = reaction.rate.falloff_coeffs
            # nsri = len(sri)
            # is_sri = True
        elif reaction.rate.type == "Lindemann":
            pass
        else:
            print(f"Unrecognized reaction rate type: {reaction.equation}")
            sys.exit(1)

    has_alpha = False
    corr_s = ""
    if not third_body and not falloff:
        pass
    elif not falloff and len(reaction.efficiencies) == 1:
        pass
    elif not falloff:
        corr_s = "alpha *"
        has_alpha = True
    else:
        corr_s = "Corr *"
        has_alpha = True

    rea_dict = OrderedDict()
    pro_dict = OrderedDict()
    all_dict = OrderedDict()
    all_wqss_dict = OrderedDict()
    sum_nuk = 0
    all_reactants = copy.deepcopy(reaction.reactants)
    all_products = copy.deepcopy(reaction.products)
    if hasattr(reaction, "efficiencies"):
        if len(reaction.efficiencies) == 1:
            all_reactants = dict(
                sum(
                    (
                        Counter(x)
                        for x in [all_reactants, reaction.efficiencies]
                    ),
                    Counter(),
                )
            )
            all_products = dict(
                sum(
                    (
                        Counter(x)
                        for x in [all_products, reaction.efficiencies]
                    ),
                    Counter(),
                )
            )

    # Build rea_dict containing reaction species
    for symbol, coefficient in all_reactants.items():
        k = species_info.ordered_idx_map[symbol]
        sum_nuk -= coefficient
        if k in rea_dict:
            coe_old = rea_dict[k][1]
            rea_dict[k] = (symbol, coefficient + coe_old)
        else:
            rea_dict[k] = (symbol, coefficient)
    # Build pro_dict containing product species
    for symbol, coefficient in all_products.items():
        k = species_info.ordered_idx_map[symbol]
        sum_nuk += coefficient
        if k in pro_dict:
            coe_old = pro_dict[k][1]
            pro_dict[k] = (symbol, coefficient + coe_old)
        else:
            pro_dict[k] = (symbol, coefficient)
    # Build the dict with species and coefficients
    for k in range(n_species):
        if k in rea_dict and k in pro_dict:
            sr, nur = rea_dict[k]
            sp, nup = pro_dict[k]
            all_dict[k] = (sr, nup - nur)
        elif k in rea_dict:
            sr, nur = rea_dict[k]
            all_dict[k] = (sr, -nur)
        elif k in pro_dict:
            sp, nup = pro_dict[k]
            all_dict[k] = (sp, nup)

    # Build the dict including qss species
    for k in range(len(species_info.all_species_list)):
        if k in rea_dict and k in pro_dict:
            sr, nur = rea_dict[k]
            sp, nup = pro_dict[k]
            all_wqss_dict[k] = (sr, nup - nur)
        elif k in rea_dict:
            sr, nur = rea_dict[k]
            all_wqss_dict[k] = (sr, -nur)
        elif k in pro_dict:
            sp, nup = pro_dict[k]
            all_wqss_dict[k] = (sp, nup)

    sorted_reactants = sorted(rea_dict.values())
    sorted_products = sorted(pro_dict.values())

    if not reaction.reversible:
        if falloff or has_alpha:
            print(
                "FIXME: irreversible reaction in _ajac_reaction may not work"
            )
            cw.writer(
                fstream,
                cw.comment(
                    "FIXME: irreversible reaction in _ajac_reaction may not work",
                ),
            )
        for k in range(n_species):
            if k in sorted_reactants and k in sorted_products:
                print(
                    "FIXME: irreversible reaction in _ajac_reaction may not work"
                )
                cw.writer(
                    fstream,
                    cw.comment(
                        "FIXME: irreversible reaction in _ajac_reaction may not work",
                    ),
                )

    if has_alpha:
        cw.writer(fstream, cw.comment("3-body correction factor"))
        enhancement_d = cu.enhancement_d(
            # mechanism, species_info, reaction, syms
            mechanism,
            species_info,
            reaction,
            syms=None,
        )
        cw.writer(fstream, "alpha = %s;" % enhancement_d)

    # forward
    cw.writer(fstream, cw.comment("forward"))
    cw.writer(
        fstream,
        "phi_f = %s;"
        % cu.qss_sorted_phase_space(
            mechanism, species_info, reaction, reaction.reactants
        ),
    )
    #
    cw.writer(fstream, "k_f = %.15g" % (pef.m))
    cw.writer(
        fstream,
        "            * exp(%.15g * tc[0] - %.15g * (%.15g) * invT);"
        % (beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, ae.m),
    )
    #
    if remove_forward:
        cw.writer(fstream, cw.comment("Remove forward reaction"))
        # DLNKFDT CHECK
        cw.writer(
            fstream,
            "dlnkfdT = %.15g * invT + %.15g * (%.15g) * invT2;"
            % (beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, ae.m),
        )
        cw.writer(fstream, cw.comment("dlnkfdT = 0.0;"))
    else:
        cw.writer(
            fstream,
            "dlnkfdT = %.15g * invT + %.15g * (%.15g) * invT2;"
            % (beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, ae.m),
        )

    if falloff:
        cw.writer(fstream, cw.comment("pressure-fall-off"))
        cw.writer(
            fstream,
            "k_0 = %.15g * exp(%.15g * tc[0] - %.15g * (%.15g) * invT);"
            % (
                low_pef.m * 10 ** (3**dim),
                low_beta,
                (1.0 / cc.Rc / cc.ureg.kelvin).m,
                low_ae.m,
            ),
        )
        cw.writer(fstream, "Pr = 1e-%d * alpha / k_f * k_0;" % (dim * 6))
        cw.writer(fstream, "fPr = Pr / (1.0+Pr);")
        cw.writer(
            fstream,
            "dlnk0dT = %.15g * invT + %.15g * (%.15g) * invT2;"
            % (low_beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, low_ae.m),
        )
        cw.writer(fstream, "dlogPrdT = log10e*(dlnk0dT - dlnkfdT);")
        cw.writer(fstream, "dlogfPrdT = dlogPrdT / (1.0+Pr);")
        #
        if is_sri:
            cw.writer(fstream, cw.comment("SRI form"))
            print("FIXME: sri not supported in _ajac_reaction yet")
            sys.exit(1)
        elif is_troe:
            cw.writer(fstream, cw.comment("Troe form"))
            troe = reaction.rate.falloff_coeffs
            ntroe = len(troe)
            cw.writer(fstream, "logPr = log10(Pr);")
            if abs(troe[1]) > 1.0e-100:
                if troe[0] < 0:
                    cw.writer(
                        fstream,
                        "Fcent1 = (1.+%.15g)*exp(-T/%.15g);"
                        % (-troe[0], troe[1]),
                    )
                else:
                    cw.writer(
                        fstream,
                        "Fcent1 = (1.-%.15g)*exp(-T/%.15g);"
                        % (troe[0], troe[1]),
                    )
            else:
                cw.writer(fstream, "Fcent1 = 0.;")
            if abs(troe[2]) > 1.0e-100:
                cw.writer(
                    fstream,
                    "Fcent2 = %.15g * exp(-T/%.15g);" % (troe[0], troe[2]),
                )
            else:
                cw.writer(fstream, "Fcent2 = 0.;")
            if ntroe == 4:
                if troe[3] < 0:
                    cw.writer(
                        fstream, "Fcent3 = exp(%.15g * invT);" % -troe[3]
                    )
                else:
                    cw.writer(
                        fstream, "Fcent3 = exp(-%.15g * invT);" % troe[3]
                    )
            else:
                cw.writer(fstream, "Fcent3 = 0.;")
            cw.writer(fstream, "Fcent = Fcent1 + Fcent2 + Fcent3;")
            cw.writer(fstream, "logFcent = log10(Fcent);")
            cw.writer(fstream, "troe_c = -.4 - .67 * logFcent;")
            cw.writer(fstream, "troe_n = .75 - 1.27 * logFcent;")
            cw.writer(
                fstream, "troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));"
            )
            cw.writer(fstream, "troePr = (troe_c + logPr) * troePr_den;")
            cw.writer(fstream, "troe = 1.0 / (1.0 + troePr*troePr);")
            cw.writer(fstream, "F = pow(10.0, logFcent * troe);")

            cw.writer(fstream, "dlogFcentdT = log10e/Fcent*( ")
            if abs(troe[1]) > 1.0e-100:
                cw.writer(fstream, "    -Fcent1/%.15g" % troe[1])
            if abs(troe[2]) > 1.0e-100:
                cw.writer(fstream, "    -Fcent2/%.15g" % troe[2])
            if ntroe == 4:
                if abs(troe[3]) > 1.0e-100:
                    cw.writer(fstream, "    + Fcent3*%.15g*invT2" % troe[3])
            cw.writer(fstream, ");")

            cw.writer(
                fstream,
                "dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr *"
                " troePr_den;",
            )
            cw.writer(
                fstream, "dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;"
            )
            cw.writer(fstream, "dlogFdn = dlogFdcn_fac * troePr;")
            cw.writer(fstream, "dlogFdlogPr = dlogFdc;")
            cw.writer(
                fstream,
                "dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) +"
                " dlogFdlogPr * dlogPrdT;",
            )
        else:
            cw.writer(fstream, cw.comment("Lindemann form"))
            cw.writer(fstream, "F = 1.0;")
            cw.writer(fstream, "dlogFdlogPr = 0.0;")
            cw.writer(fstream, "dlogFdT = 0.0;")

    # reverse
    if not reaction.reversible:
        cw.writer(fstream, cw.comment("rate of progress"))
        if (not has_alpha) and (not falloff):
            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(fstream, cw.comment("q = k_f*phi_f;"))
                cw.writer(fstream, "q = 0.0;")
            else:
                cw.writer(fstream, "q = k_f*phi_f;")
        else:
            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(fstream, cw.comment("q_nocor = k_f*phi_f;"))
                cw.writer(fstream, "q_nocor = 0.0;")
            else:
                cw.writer(fstream, "q_nocor = k_f*phi_f;")

            if falloff:
                cw.writer(fstream, "Corr = fPr * F;")
                cw.writer(fstream, "q = Corr * q_nocor;")
            else:
                cw.writer(fstream, "q = alpha * q_nocor;")

        if falloff:
            cw.writer(fstream, "dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);")
            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(
                    fstream,
                    cw.comment(
                        "dqdT = %sdlnkfdT*k_f*phi_f + dlnCorrdT*q;" % corr_s,
                    ),
                )
                cw.writer(fstream, "dqdT = dlnCorrdT*q;")
            else:
                cw.writer(
                    fstream,
                    "dqdT = %sdlnkfdT*k_f*phi_f + dlnCorrdT*q;" % corr_s,
                )
        else:
            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(
                    fstream, cw.comment("dqdT = %sdlnkfdT*k_f*phi_f;" % corr_s)
                )
                cw.writer(fstream, "dqdT = 0;")
            else:
                cw.writer(fstream, "dqdT = %sdlnkfdT*k_f*phi_f;" % corr_s)
    else:
        cw.writer(fstream, cw.comment("reverse"))
        cw.writer(
            fstream,
            "phi_r = %s;"
            % cu.qss_sorted_phase_space(
                mechanism, species_info, reaction, reaction.products
            ),
        )
        cw.writer(
            fstream,
            "Kc = %s;" % cu.sorted_kc(mechanism, species_info, reaction),
        )
        cw.writer(fstream, "k_r = k_f / Kc;")

        dlnkcdt_s = "invT * ("
        terms = []
        for symbol, coefficient in sorted(
            sorted_reactants, key=lambda v: species_info.dict_species[v[0]]
        ):
            k = species_info.ordered_idx_map[symbol]
            if symbol not in species_info.qssa_species_list:
                if coefficient == 1.0:
                    terms.append("h_RT[%d]" % (k))
                else:
                    terms.append("%f*h_RT[%d]" % (coefficient, k))
            else:
                if coefficient == 1.0:
                    terms.append("h_RT_qss[%d]" % (k - n_species))
                else:
                    terms.append(
                        "%f*h_RT_qss[%d]" % (coefficient, k - n_species)
                    )
        dlnkcdt_s += "-(" + " + ".join(terms) + ")"
        terms = []
        for symbol, coefficient in sorted(
            sorted_products, key=lambda v: species_info.dict_species[v[0]]
        ):
            k = species_info.ordered_idx_map[symbol]
            if symbol not in species_info.qssa_species_list:
                if coefficient == 1.0:
                    terms.append("h_RT[%d]" % (k))
                else:
                    terms.append("%f*h_RT[%d]" % (coefficient, k))
            else:
                if coefficient == 1.0:
                    terms.append("h_RT_qss[%d]" % (k - n_species))
                else:
                    terms.append(
                        "%f*h_RT_qss[%d]" % (coefficient, k - n_species)
                    )
        dlnkcdt_s += " + (" + " + ".join(terms) + ")"
        if sum_nuk > 0:
            dlnkcdt_s += " - %f" % sum_nuk
        elif sum_nuk < 0:
            dlnkcdt_s += " + %f" % (-sum_nuk)
        dlnkcdt_s += ")"
        cw.writer(fstream, "dlnKcdT = %s;" % dlnkcdt_s)

        cw.writer(fstream, "dkrdT = (dlnkfdT - dlnKcdT)*k_r;")

        cw.writer(fstream, cw.comment("rate of progress"))
        if (not has_alpha) and (not falloff):
            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(fstream, cw.comment("q = k_f*phi_f - k_r*phi_r;"))
                cw.writer(fstream, "q = - k_r*phi_r;")
            else:
                cw.writer(fstream, "q = k_f*phi_f - k_r*phi_r;")
        else:
            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(
                    fstream, cw.comment("q_nocor = k_f*phi_f - k_r*phi_r;")
                )
                cw.writer(fstream, "q_nocor = - k_r*phi_r;")
            else:
                cw.writer(fstream, "q_nocor = k_f*phi_f - k_r*phi_r;")
            if falloff:
                cw.writer(fstream, "Corr = fPr * F;")
                cw.writer(fstream, "q = Corr * q_nocor;")
            else:
                cw.writer(fstream, "q = alpha * q_nocor;")

        if falloff:
            cw.writer(fstream, "dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);")
            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(
                    fstream,
                    cw.comment(
                        "dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r) +"
                        " dlnCorrdT*q;" % corr_s,
                    ),
                )
                cw.writer(
                    fstream, "dqdT = %s(- dkrdT*phi_r) + dlnCorrdT*q;" % corr_s
                )
            else:
                cw.writer(
                    fstream,
                    "dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;"
                    % corr_s,
                )
        else:
            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(
                    fstream,
                    cw.comment(
                        "dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r);" % corr_s,
                    ),
                )
                cw.writer(fstream, "dqdT = %s( - dkrdT*phi_r);" % corr_s)
            else:
                cw.writer(
                    fstream,
                    "dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r);" % corr_s,
                )

    cw.writer(fstream, cw.comment("update wdot"))
    for k in sorted(all_dict.keys()):
        s, nu = all_dict[k]
        if nu == 1:
            cw.writer(fstream, "wdot[%d] += q;" % (k) + cw.comment("%s" % (s)))
        elif nu == -1:
            cw.writer(fstream, "wdot[%d] -= q;" % (k) + cw.comment("%s" % (s)))
        elif nu > 0:
            cw.writer(
                fstream,
                "wdot[%d] += %.15g * q;" % (k, nu) + cw.comment("%s" % (s)),
            )
        elif nu < 0:
            cw.writer(
                fstream,
                "wdot[%d] -= %.15g * q;" % (k, -nu) + cw.comment("%s" % (s)),
            )

    if falloff:
        cw.writer(fstream, cw.comment("for convenience"))
        cw.writer(fstream, "k_f *= Corr;")
        if reaction.reversible:
            cw.writer(fstream, "k_r *= Corr;")
    elif has_alpha:
        cw.writer(fstream, cw.comment("for convenience"))
        cw.writer(fstream, "k_f *= alpha;")
        if reaction.reversible:
            cw.writer(fstream, "k_r *= alpha;")
        else:
            cw.writer(fstream, "k_r = 0.0;")

    if falloff:
        if precond:
            cw.writer(fstream, "dcdc_fac = 0.0;")
        else:
            cw.writer(
                fstream, "dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);"
            )
    # elif has_alpha:
    #    cw.writer(fstream,'dcdc_fac = q_nocor;')

    if has_alpha or falloff:

        if not precond:
            cw.writer(fstream, "if (consP) {")

            for k in range(n_species):
                dqdc_s = denhancement_d(
                    mechanism, species_info, reaction, k, True
                )
                if dqdc_s != "0":
                    if falloff:
                        if dqdc_s == "1":
                            dqdc_s = "dcdc_fac"
                        else:
                            dqdc_s += "*dcdc_fac"
                    elif has_alpha:
                        if dqdc_s == "1":
                            dqdc_s = "q_nocor"
                        else:
                            dqdc_s += "*q_nocor"

                dqdc_s = dqdc_d(
                    fstream,
                    mechanism,
                    species_info,
                    reaction,
                    sorted_reactants,
                    sorted_products,
                    rea_dict,
                    pro_dict,
                    dqdc_s,
                    k,
                    remove_forward,
                    syms,
                )
                if dqdc_s:
                    symb_k = species_info.all_species_list[k]
                    cw.writer(fstream, cw.comment("d()/d[%s]" % symb_k))
                    cw.writer(fstream, "dqdci = %s;" % (dqdc_s))
                    #
                    for m in sorted(all_dict.keys()):
                        if all_dict[m][1] != 0:
                            s1 = "J[%d] += %.15g * dqdci;" % (
                                k * (n_species + 1) + m,
                                all_dict[m][1],
                            )
                            s1 = s1.replace("+= 1 *", "+=").replace(
                                "+= -1 *", "-="
                            )
                            s2 = cw.comment(
                                "dwdot[%s]/d[%s]"
                                % (
                                    all_dict[m][0],
                                    symb_k,
                                )
                            )
                            cw.writer(fstream, s1.ljust(30) + s2)

            cw.writer(fstream, "}")
            cw.writer(fstream, "else {")

        for k in range(n_species):
            # for k in range(len(species_info.all_species_list)):
            dqdc_s = denhancement_d(
                mechanism, species_info, reaction, k, False
            )
            if dqdc_s != "0":
                if falloff:
                    if dqdc_s == "1":
                        dqdc_s = "dcdc_fac"
                    else:
                        dqdc_s += "*dcdc_fac"
                elif has_alpha:
                    if dqdc_s == "1":
                        dqdc_s = "q_nocor"
                    else:
                        dqdc_s += "*q_nocor"

            dqdc_s = dqdc_d(
                fstream,
                mechanism,
                species_info,
                reaction,
                sorted_reactants,
                sorted_products,
                rea_dict,
                pro_dict,
                dqdc_s,
                k,
                remove_forward,
                syms,
            )
            if dqdc_s:
                cw.writer(fstream, "dqdc[%d] = %s;" % (k, dqdc_s))
            elif precond:
                cw.writer(fstream, "dqdc[%d] = 0.0;" % k)

        cw.writer(fstream, "for (int k=0; k<%d; k++) {" % n_species)
        for m in sorted(all_dict.keys()):
            if all_dict[m][1] != 0:
                s1 = "J[%d*k+%d] += %.15g * dqdc[k];" % (
                    (n_species + 1),
                    m,
                    all_dict[m][1],
                )
                s1 = s1.replace("+= 1 *", "+=").replace("+= -1 *", "-=")
                cw.writer(fstream, s1)
        cw.writer(fstream, "}")

        if not precond:
            cw.writer(fstream, "}")

        for m in sorted(all_dict.keys()):
            if all_dict[m][1] != 0:
                s1 = "J[%d] += %.15g * dqdT;" % (
                    n_species * (n_species + 1) + m,
                    all_dict[m][1],
                ) + cw.comment("dwdot[%s]/dT" % (all_dict[m][0]))
                s1 = s1.replace("+= 1 *", "+=").replace("+= -1 *", "-=")
                cw.writer(fstream, s1)

    else:

        for k in range(n_species):
            dqdc_s = dqdc_d(
                fstream,
                mechanism,
                species_info,
                reaction,
                sorted_reactants,
                sorted_products,
                rea_dict,
                pro_dict,
                "",
                k,
                remove_forward,
                syms,
            )

            if dqdc_s:
                cw.writer(
                    fstream, cw.comment("d()/d[%s]" % all_wqss_dict[k][0])
                )
                cw.writer(fstream, "dqdci = %s;" % (dqdc_s))
                if reaction.reversible or k in rea_dict:
                    for m in sorted(all_dict.keys()):
                        if all_dict[m][1] != 0:
                            s1 = "J[%d] += %.15g * dqdci;" % (
                                k * (n_species + 1) + m,
                                all_wqss_dict[m][1],
                            )
                            s1 = s1.replace("+= 1 *", "+=").replace(
                                "+= -1 *", "-="
                            )
                            s2 = cw.comment(
                                "dwdot[%s]/d[%s]"
                                % (
                                    all_wqss_dict[m][0],
                                    all_wqss_dict[k][0],
                                )
                            )
                            cw.writer(fstream, s1.ljust(30) + s2)
        cw.writer(fstream, cw.comment("d()/dT"))
        for m in sorted(all_dict.keys()):
            if all_dict[m][1] != 0:
                s1 = "J[%d] += %.15g * dqdT;" % (
                    n_species * (n_species + 1) + m,
                    all_dict[m][1],
                )
                s1 = (
                    s1.replace("+= 1 *", "+=")
                    .replace("+= -1 *", "-=")
                    .replace("+= -1 *", "-=")
                )
                s2 = cw.comment("dwdot[%s]/dT" % (all_dict[m][0]))
                cw.writer(fstream, s1.ljust(30) + s2)


def dqdc_d(
    fstream,
    mechanism,
    species_info,
    reaction,
    sorted_reactants,
    sorted_products,
    rea_dict,
    pro_dict,
    dqdc_s,
    k,
    remove_forward,
    syms,
):
    """Write dqdc."""
    if dqdc_s == "0":
        dqdc_s = ""
    if k in sorted(rea_dict.keys()):
        dps = dphase_space(
            mechanism,
            species_info,
            sorted_reactants,
            rea_dict[k][0],
            syms,
        )
        if dps == "1.0":
            dps_s = ""
        else:
            dps_s = "*" + dps
        if remove_forward:
            cw.writer(fstream, cw.comment("Remove forward reaction"))
            dqdc_s += ""
        else:
            dqdc_s += " + k_f%s" % dps_s
    if reaction.reversible:
        if k in sorted(pro_dict.keys()):
            dps = dphase_space(
                mechanism,
                species_info,
                sorted_products,
                pro_dict[k][0],
                syms,
            )
            if dps == "1.0":
                dps_s = ""
            else:
                dps_s = "*" + dps
            dqdc_s += " - k_r%s" % dps_s
    return dqdc_s


def denhancement_d(mechanism, species_info, reaction, kid, cons_p):
    """Get enhancement gradient."""
    third_body = reaction.reaction_type == "three-body"
    falloff = reaction.reaction_type == "falloff"
    if not third_body and not falloff:
        print("denhancement_d called for a reaction without a third body")
        sys.exit(1)

    if not hasattr(reaction, "efficiencies"):
        print("FIXME")
        sys.exit(1)
        species, coefficient = third_body
        if species == "<mixture>":
            if cons_p:
                return "0"
            else:
                return "1"
        elif species_info.ordered_idx_map[species] == kid:
            return "1"
        else:
            return "0"
    else:
        efficiencies = reaction.efficiencies
        if cons_p:
            for _, (symbol, efficiency) in enumerate(efficiencies.items()):
                if species_info.ordered_idx_map[symbol] == kid:
                    return "(%.15g - 1)" % (efficiency)
            return "0"
        else:
            for _, (symbol, efficiency) in enumerate(efficiencies.items()):
                if species_info.ordered_idx_map[symbol] == kid:
                    return "%.15g" % (efficiency)
            return "1"


def dphase_space(mechanism, species_info, reagents, r, syms):
    """Get string for phase space gradient."""
    reagents = {x[0]: x[1] for x in reagents}
    phi = []
    sorted_reagents = sorted(
        reagents.keys(), key=lambda v: species_info.dict_species[v]
    )
    for symbol in sorted_reagents:
        coefficient = reagents[symbol]
        if symbol not in species_info.qssa_species_list:
            if symbol == r:
                if coefficient > 1:
                    phi += ["%f" % coefficient]
                    if (coefficient - 1) == 1.0:
                        conc = "sc[%d]" % (
                            species_info.ordered_idx_map[symbol]
                        )
                    else:
                        conc = "pow(sc[%d],%f)" % (
                            species_info.ordered_idx_map[symbol],
                            (coefficient - 1),
                        )
                    phi += [conc]
            else:
                if coefficient == 1.0:
                    conc = "sc[%d]" % (species_info.ordered_idx_map[symbol])
                else:
                    conc = "pow(sc[%d], %f)" % (
                        species_info.ordered_idx_map[symbol],
                        coefficient,
                    )
                phi += [conc]
        # Symbol is in qssa_species_list
        else:
            if symbol == r:
                if coefficient > 1:
                    phi += ["%f" % coefficient]
                    if (coefficient - 1) == 1.0:
                        conc = "sc_qss[%d]" % (
                            species_info.ordered_idx_map[symbol]
                            - species_info.n_species
                        )
                    else:
                        conc = "pow(sc_qss[%d],%f)" % (
                            species_info.ordered_idx_map[symbol]
                            - species_info.n_species,
                            (coefficient - 1),
                        )
                    phi += [conc]
            else:
                if coefficient == 1.0:
                    conc = "sc_qss[%d]" % (
                        species_info.ordered_idx_map[symbol]
                        - species_info.n_species
                    )
                else:
                    conc = "pow(sc_qss[%d], %f)" % (
                        species_info.ordered_idx_map[symbol]
                        - species_info.n_species,
                        coefficient,
                    )
                phi += [conc]

    if phi:
        return "*".join(phi)
    else:
        return "1.0"


def dproduction_rate(
    fstream, mechanism, species_info, reaction_info, precond=False
):
    """Write the reaction jacobian."""
    n_species = species_info.n_species

    cw.writer(fstream)
    if precond:
        cw.writer(
            fstream,
            cw.comment(
                "compute an approx to the reaction Jacobian (for preconditioning)"
            ),
        )
        cw.writer(
            fstream,
            "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void"
            " DWDOT_SIMPLIFIED(amrex::Real *  J, amrex::Real *  sc, amrex::Real * "
            " Tp, const int * HP)",
        )
    else:
        cw.writer(fstream, cw.comment("compute the reaction Jacobian"))
        cw.writer(
            fstream,
            "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *  J,"
            " amrex::Real *  sc, amrex::Real *  Tp, const int * consP)",
        )

    cw.writer(fstream, "{")
    cw.writer(fstream, "amrex::Real c[%d];" % (n_species))
    cw.writer(fstream)
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % n_species)
    cw.writer(fstream, "c[k] = 1.e6 * sc[k];")
    cw.writer(fstream, "}")

    cw.writer(fstream)
    if precond:
        cw.writer(fstream, "aJacobian_precond(J, c, *Tp, *HP);")
    else:
        cw.writer(fstream, "aJacobian(J, c, *Tp, *consP);")

    cw.writer(fstream)
    cw.writer(fstream, cw.comment("dwdot[k]/dT"))
    cw.writer(fstream, cw.comment("dTdot/d[X]"))
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % n_species)
    cw.writer(fstream, "J[%d+k] *= 1.e-6;" % (n_species * (n_species + 1)))
    cw.writer(fstream, "J[k*%d+%d] *= 1.e6;" % (n_species + 1, n_species))
    cw.writer(fstream, "}")

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
    dwdotdc=None,
    index=0,
):
    """Temporary Write jacobian term for debugging."""
    n_species = species_info.n_species
    n_reactions = mechanism.n_reactions

    cw.writer(fstream)

    # main
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
        " ajac_term_debug(amrex::Real * J, amrex::Real * sc, amrex::Real T)",
    )
    cw.writer(fstream, "{")

    # Initialize the big Jacobian array
    cw.writer(fstream, "for (int i=0; i<%d; i++) {" % (n_species + 1) ** 2)
    cw.writer(fstream, "J[i] = 0.0;")
    cw.writer(fstream, "}")
   
    cw.writer(fstream, cw.comment(f"J corresponds to index: {index}"))
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

            syms.write_array_to_cpp([dwdotdc], f"J", cw, fstream)

            cw.writer(fstream)


    # dwdotdT
    cw.writer(
        fstream,
        "amrex::Real T_pert1, T_pert2, pertT;"
    )
    cw.writer(
        fstream,
        "amrex::Real wdot_pert1[%d], wdot_pert2[%d];"
        % (
            n_species,
            n_species,
        ),
    )
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("dwdot/dT by finite difference"))
    cw.writer(fstream,"pertT = 1e-1;") 
    cw.writer(fstream,"T_pert1 = T + pertT;")
    cw.writer(fstream,"T_pert2 = T - pertT;")
    cw.writer(fstream)
    cw.writer(fstream,"productionRate(wdot_pert1, sc, T_pert1);")
    cw.writer(fstream,"productionRate(wdot_pert2, sc, T_pert2);")
    cw.writer(fstream)
    cw.writer(fstream, "for (int k = 0; k < %d ; k++) {" % n_species)
    cw.writer(
        fstream, 
        "J[%d + k] = (wdot_pert1[k] - wdot_pert2[k])/(2.0*pertT);" 
        % (
            n_species*(n_species + 1),
        )  
    )
    cw.writer(fstream, "}")

    cw.writer(fstream)

    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")

    cw.writer(fstream)
