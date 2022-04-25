import sys
from collections import OrderedDict

import ceptr.constants as cc
import ceptr.utilities as cu
import ceptr.writer as cw


def ajacPrecond(fstream, mechanism, species_info, reaction_info):
    nSpecies = species_info.nSpecies
    nReactions = reaction_info.nReactions

    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("compute an approx to the reaction Jacobian")
    )

    # cw.writer(fstream,"#ifdef COMPILE_JACOBIAN")
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void"
        " aJacobian_precond(amrex::Real *  J, amrex::Real *  sc, amrex::Real"
        " T, int HP)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "for (int i=0; i<%d; i++) {" % (nSpecies + 1) ** 2)
    cw.writer(fstream, "J[i] = 0.0;")
    cw.writer(fstream, "}")

    cw.writer(fstream)

    cw.writer(fstream, "amrex::Real wdot[%d];" % (nSpecies))
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies))
    cw.writer(fstream, "wdot[k] = 0.0;")
    cw.writer(fstream, "}")

    cw.writer(fstream)

    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T }; //"
        " temperature cache",
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
    cw.writer(fstream, "for (int k = 0; k < %d; ++k) {" % nSpecies)
    cw.writer(fstream, "mixture += sc[k];")
    cw.writer(fstream, "}")

    cw.writer(fstream)

    cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
    cw.writer(fstream, "amrex::Real g_RT[%d];" % (nSpecies))
    cw.writer(fstream, "gibbs(g_RT, tc);")
    if species_info.nQSSspecies > 0:
        cw.writer(
            fstream, "amrex::Real g_RT_qss[%d];" % (species_info.nQSSspecies)
        )
        cw.writer(fstream, "gibbs_qss(g_RT_qss, tc);")

    cw.writer(fstream)

    cw.writer(fstream, cw.comment("compute the species enthalpy"))
    cw.writer(fstream, "amrex::Real h_RT[%d];" % (nSpecies))
    cw.writer(fstream, "speciesEnthalpy(h_RT, tc);")
    if species_info.nQSSspecies > 0:
        cw.writer(
            fstream, "amrex::Real h_RT_qss[%d];" % (species_info.nQSSspecies)
        )
        cw.writer(fstream, "speciesEnthalpy_qss(h_RT_qss, tc);")

    if species_info.nQSSspecies > 0:
        cw.writer(fstream, "// Fill sc_qss here")
        cw.writer(
            fstream, "amrex::Real sc_qss[%d];" % species_info.nQSSspecies
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
    cw.writer(fstream, "amrex::Real dqdci, dcdc_fac, dqdc[%d];" % (nSpecies))
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

    for orig_idx, idx in reaction_info.idxmap.items():
        reaction = mechanism.reaction(orig_idx)
        # lt = reaction.lt
        # if lt:
        #     print("Landau-Teller reactions are not supported")
        #     sys.exit(1)

        cw.writer(
            fstream,
            cw.comment("reaction %d: %s" % (orig_idx, reaction.equation)),
        )
        falloff = reaction.reaction_type == "falloff"
        thirdBody = reaction.reaction_type == "three-body"
        if falloff:  # case 1
            cw.writer(fstream, cw.comment("a pressure-fall-off reaction"))
            ajac_reaction_precond(
                fstream,
                mechanism,
                species_info,
                reaction_info,
                reaction,
                orig_idx,
                1,
            )
        elif thirdBody:  # case 2
            cw.writer(
                fstream,
                cw.comment("a third-body and non-pressure-fall-off reaction"),
            )
            ajac_reaction_precond(
                fstream,
                mechanism,
                species_info,
                reaction_info,
                reaction,
                orig_idx,
                2,
            )
        else:  # case 3
            cw.writer(
                fstream,
                cw.comment(
                    "a non-third-body and non-pressure-fall-off reaction"
                ),
            )
            ajac_reaction_precond(
                fstream,
                mechanism,
                species_info,
                reaction_info,
                reaction,
                orig_idx,
                3,
            )
        cw.writer(fstream)

    cw.writer(
        fstream,
        "amrex::Real c_R[%d], dcRdT[%d], e_RT[%d];"
        % (nSpecies, nSpecies, nSpecies),
    )
    cw.writer(fstream, "amrex::Real * eh_RT;")
    cw.writer(fstream, "if (HP) {")

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
    cw.writer(fstream, "for (int k = 0; k < %d; ++k) {" % nSpecies)
    cw.writer(fstream, "cmix += c_R[k]*sc[k];")
    cw.writer(fstream, "dcmixdT += dcRdT[k]*sc[k];")
    cw.writer(fstream, "ehmix += eh_RT[k]*wdot[k];")
    cw.writer(
        fstream,
        "dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[%d+k];"
        % (nSpecies * (nSpecies + 1)),
    )
    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "amrex::Real cmixinv = 1.0/cmix;")
    cw.writer(fstream, "amrex::Real tmp1 = ehmix*cmixinv;")
    cw.writer(fstream, "amrex::Real tmp3 = cmixinv*T;")
    cw.writer(fstream, "amrex::Real tmp2 = tmp1*tmp3;")
    cw.writer(fstream, "amrex::Real dehmixdc;")

    cw.writer(fstream, "// dTdot/d[X]")
    cw.writer(fstream, "for (int k = 0; k < %d; ++k) {" % nSpecies)
    cw.writer(fstream, "dehmixdc = 0.0;")
    cw.writer(fstream, "for (int m = 0; m < %d; ++m) {" % nSpecies)
    cw.writer(fstream, "dehmixdc += eh_RT[m]*J[k*%s+m];" % (nSpecies + 1))
    cw.writer(fstream, "}")
    cw.writer(
        fstream,
        "J[k*%d+%d] = tmp2*c_R[k] - tmp3*dehmixdc;" % (nSpecies + 1, nSpecies),
    )
    cw.writer(fstream, "}")

    cw.writer(fstream, "// dTdot/dT")
    cw.writer(
        fstream,
        "J[%d] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;"
        % (nSpecies * (nSpecies + 1) + nSpecies),
    )

    cw.writer(fstream, "}")


def ajac_reaction_precond(
    fstream, mechanism, species_info, reaction_info, reaction, orig_idx, rcase
):

    nSpecies = species_info.nSpecies
    removeForward = cu.isRemoveForward(reaction_info, orig_idx)
    thirdBody = reaction.reaction_type == "three-body"
    if rcase == 1:  # pressure-dependent reaction
        isPD = True
        # FIXME is this right?
        has_alpha = True
        cw.writer(fstream, "// also 3-body")
        # if thirdBody:
        #     has_alpha = True
        #     cw.writer(fstream,"// also 3-body")
        # else:
        #     has_alpha = False
        #     cw.writer(fstream,"// non 3-body")
        #     # print(
        #     #     "FIXME: pressure dependent non-3-body reaction in _ajac_reaction"
        #     # )
        #     # sys.exit(1)
    elif rcase == 2:  # third-body and non-pressure-dependent reaction
        isPD = False
        has_alpha = True
    elif rcase == 3:  # simple non-third and non-pressure-dependent reaction
        isPD = False
        has_alpha = False
    else:
        print("_ajac_reaction: wrong case ", rcase)
        exit(1)

    rea_dict = {}
    pro_dict = {}
    all_dict = {}
    sumNuk = 0
    for symbol, coefficient in reaction.reactants.items():
        k = species_info.ordered_idx_map[symbol]
        sumNuk -= coefficient
        if k in rea_dict:
            coe_old = rea_dict[k][1]
            rea_dict[k] = (symbol, coefficient + coe_old)
        else:
            rea_dict[k] = (symbol, coefficient)
    for symbol, coefficient in reaction.products.items():
        k = species_info.ordered_idx_map[symbol]
        sumNuk += coefficient
        if k in pro_dict:
            coe_old = pro_dict[k][1]
            pro_dict[k] = (symbol, coefficient + coe_old)
        else:
            pro_dict[k] = (symbol, coefficient)
    for k in range(nSpecies):
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

    sorted_reactants = sorted(rea_dict.values())
    sorted_products = sorted(pro_dict.values())

    if not reaction.reversible:
        if isPD or has_alpha:
            print(
                "FIXME: irreversible reaction in _ajac_reaction may not work"
            )
            cw.writer(
                fstream,
                "// FIXME: irreversible reaction in _ajac_reaction may not"
                " work",
            )
        for k in range(nSpecies):
            if k in sorted_reactants and k in sorted_products:
                print(
                    "FIXME: irreversible reaction in _ajac_reaction may not"
                    " work"
                )
                cw.writer(
                    fstream,
                    "// FIXME: irreversible reaction in _ajac_reaction may not"
                    " work",
                )

    dim = cu.phaseSpaceUnits(reaction.reactants)
    aeuc = cu.activationEnergyUnits()
    is_sri = False
    is_troe = False
    if isPD:
        Corr_s = "Corr *"
        uc = cu.prefactorUnits(
            cc.ureg("mole/cm**3"), 1 - dim
        )  # Case 2 !PD, TB
        ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), 1 - dim)
        A = (reaction.high_rate.pre_exponential_factor * ctuc).to(
            uc.to_root_units()
        )
        beta = reaction.high_rate.temperature_exponent
        E = (
            reaction.high_rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
        ).to(aeuc)

        low_A = (reaction.low_rate.pre_exponential_factor * ctuc).to(
            uc.to_root_units()
        )
        low_beta = reaction.low_rate.temperature_exponent
        low_E = (
            reaction.low_rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
        ).to(aeuc)
        if reaction.rate.type == "Troe":
            troe = reaction.rate.falloff_coeffs
            ntroe = len(troe)
            is_troe = True
        elif reaction.rate.type == "Sri":
            sri = reaction.rate.falloff_coeffs
            nsri = len(sri)
            is_sri = True
        else:
            print("Unrecognized reaction rate type")
            sys.exit(1)
    elif has_alpha:
        Corr_s = "alpha * "
        uc = cu.prefactorUnits(cc.ureg("mole/cm**3"), -dim)  # Case 2 !PD, TB
        ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), -dim)
        A = (reaction.rate.pre_exponential_factor * ctuc).to(
            uc.to_root_units()
        )
        beta = reaction.rate.temperature_exponent
        E = (
            reaction.rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
        ).to(aeuc)
    else:
        uc = cu.prefactorUnits(
            cc.ureg("mole/cm**3"), 1 - dim
        )  # Case 3 !PD, !TB
        ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), 1 - dim)
        A = (reaction.rate.pre_exponential_factor * ctuc).to(
            uc.to_root_units()
        )
        beta = reaction.rate.temperature_exponent
        E = (
            reaction.rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
        ).to(aeuc)
        Corr_s = ""

    if has_alpha:
        cw.writer(fstream, "// 3-body correction factor")
        cw.writer(
            fstream,
            "alpha = %s;" % enhancement_d(mechanism, species_info, reaction),
        )

    # forward
    cw.writer(fstream, "// forward")
    cw.writer(
        fstream,
        "phi_f = %s;"
        % cu.QSSsortedPhaseSpace(mechanism, species_info, reaction.reactants),
    )
    #
    cw.writer(fstream, "k_f = %.15g" % (A.m))
    cw.writer(
        fstream,
        "            * exp(%.15g * tc[0] - %.15g * (%.15g) * invT);"
        % (beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, E.m),
    )
    #
    if removeForward:
        cw.writer(fstream, "// Remove forward reaction")
        # DLNKFDT CHECK
        cw.writer(
            fstream,
            "dlnkfdT = %.15g * invT + %.15g *  (%.15g)  * invT2;"
            % (beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, E.m),
        )
        cw.writer(fstream, "//dlnkfdT = 0.0;")
    else:
        cw.writer(
            fstream,
            "dlnkfdT = %.15g * invT + %.15g *  (%.15g)  * invT2;"
            % (beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, E.m),
        )

    if isPD:
        cw.writer(fstream, "// pressure-fall-off")
        cw.writer(
            fstream,
            "k_0 = %.15g * exp(%.15g * tc[0] - %.15g * (%.15g) * invT);"
            % (
                low_A.m * 10 ** (3**dim),
                low_beta,
                (1.0 / cc.Rc / cc.ureg.kelvin).m,
                low_E.m,
            ),
        )
        cw.writer(fstream, "Pr = 1e-%d * alpha / k_f * k_0;" % (dim * 6))
        cw.writer(fstream, "fPr = Pr / (1.0+Pr);")
        cw.writer(
            fstream,
            "dlnk0dT = %.15g * invT + %.15g * (%.15g) * invT2;"
            % (low_beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, low_E.m),
        )
        cw.writer(fstream, "dlogPrdT = log10e*(dlnk0dT - dlnkfdT);")
        cw.writer(fstream, "dlogfPrdT = dlogPrdT / (1.0+Pr);")
        #
        if is_sri:
            cw.writer(fstream, "// SRI form")
            print("FIXME: sri not supported in _ajac_reaction yet")
            sys.exit(1)
        elif is_troe:
            cw.writer(fstream, "// Troe form")
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
            cw.writer(fstream, "// Lindemann form")
            cw.writer(fstream, "F = 1.0;")
            cw.writer(fstream, "dlogFdlogPr = 0.0;")
            cw.writer(fstream, "dlogFdT = 0.0;")

    # reverse
    if not reaction.reversible:
        cw.writer(fstream, "// rate of progress")
        if (not has_alpha) and (not isPD):
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(fstream, "//q = k_f*phi_f;")
                cw.writer(fstream, "q = 0.0;")
            else:
                cw.writer(fstream, "q = k_f*phi_f;")
        else:
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(fstream, "//q_nocor = k_f*phi_f;")
                cw.writer(fstream, "q_nocor = 0.0;")
            else:
                cw.writer(fstream, "q_nocor = k_f*phi_f;")
            if isPD:
                cw.writer(fstream, "Corr = fPr * F;")
                cw.writer(fstream, "q = Corr * q_nocor;")
            else:
                cw.writer(fstream, "q = alpha * q_nocor;")

        if isPD:
            cw.writer(fstream, "dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);")
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(
                    fstream,
                    "//dqdT = %sdlnkfdT*k_f*phi_f + dlnCorrdT*q;" % Corr_s,
                )
                cw.writer(fstream, "dqdT = dlnCorrdT*q;")
            else:
                cw.writer(
                    fstream,
                    "dqdT = %sdlnkfdT*k_f*phi_f + dlnCorrdT*q;" % Corr_s,
                )
        else:
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(fstream, "//dqdT = %sdlnkfdT*k_f*phi_f;" % Corr_s)
                cw.writer(fstream, "dqdT = 0;")
            else:
                cw.writer(fstream, "dqdT = %sdlnkfdT*k_f*phi_f;" % Corr_s)
    else:
        cw.writer(fstream, "// reverse")
        cw.writer(
            fstream,
            "phi_r = %s;"
            % cu.QSSsortedPhaseSpace(
                mechanism, species_info, reaction.products
            ),
        )
        cw.writer(
            fstream,
            "Kc = %s;" % cu.sortedKc(mechanism, species_info, reaction),
        )
        cw.writer(fstream, "k_r = k_f / Kc;")

        dlnKcdT_s = "invT * ("
        terms = []
        for symbol, coefficient in sorted(
            sorted_reactants,
            key=lambda x: next(
                (y for y in species_info.all_species if y.name == x[0]),
                None,
            ).idx,
        ):
            k = species_info.ordered_idx_map[symbol]
            if symbol not in species_info.qss_species_list:
                if coefficient == 1.0:
                    terms.append("h_RT[%d]" % (k))
                else:
                    terms.append("%f*h_RT[%d]" % (coefficient, k))
            else:
                if coefficient == 1.0:
                    terms.append("h_RT_qss[%d]" % (k - nSpecies))
                else:
                    terms.append(
                        "%f*h_RT_qss[%d]" % (coefficient, k - nSpecies)
                    )
        dlnKcdT_s += "-(" + " + ".join(terms) + ")"
        terms = []
        for symbol, coefficient in sorted(
            sorted_products,
            key=lambda x: next(
                (y for y in species_info.all_species if y.name == x[0]),
                None,
            ).idx,
        ):
            k = species_info.ordered_idx_map[symbol]
            if symbol not in species_info.qss_species_list:
                if coefficient == 1.0:
                    terms.append("h_RT[%d]" % (k))
                else:
                    terms.append("%f*h_RT[%d]" % (coefficient, k))
            else:
                if coefficient == 1.0:
                    terms.append("h_RT_qss[%d]" % (k - nSpecies))
                else:
                    terms.append(
                        "%f*h_RT_qss[%d]" % (coefficient, k - nSpecies)
                    )
        dlnKcdT_s += " + (" + " + ".join(terms) + ")"
        if sumNuk > 0:
            dlnKcdT_s += " - %f" % sumNuk
        elif sumNuk < 0:
            dlnKcdT_s += " + %f" % (-sumNuk)
        dlnKcdT_s += ")"
        cw.writer(fstream, "dlnKcdT = %s;" % dlnKcdT_s)

        cw.writer(fstream, "dkrdT = (dlnkfdT - dlnKcdT)*k_r;")

        cw.writer(fstream, "// rate of progress")
        if (not has_alpha) and (not isPD):
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(fstream, "//q = k_f*phi_f - k_r*phi_r;")
                cw.writer(fstream, "q = - k_r*phi_r;")
            else:
                cw.writer(fstream, "q = k_f*phi_f - k_r*phi_r;")
        else:
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(fstream, "//q = k_f*phi_f - k_r*phi_r;")
                cw.writer(fstream, "q_nocor = - k_r*phi_r;")
            else:
                cw.writer(fstream, "q_nocor = k_f*phi_f - k_r*phi_r;")
            if isPD:
                cw.writer(fstream, "Corr = fPr * F;")
                cw.writer(fstream, "q = Corr * q_nocor;")
            else:
                cw.writer(fstream, "q = alpha * q_nocor;")

        if isPD:
            cw.writer(fstream, "dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);")
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(
                    fstream,
                    "//dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r) +"
                    " dlnCorrdT*q;" % Corr_s,
                )
                cw.writer(
                    fstream,
                    "dqdT = %s( - dkrdT*phi_r) + dlnCorrdT*q;" % Corr_s,
                )
            else:
                cw.writer(
                    fstream,
                    "dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;"
                    % Corr_s,
                )
        else:
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(
                    fstream,
                    "//dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r);" % Corr_s,
                )
                cw.writer(fstream, "dqdT = %s( - dkrdT*phi_r);" % Corr_s)
            else:
                cw.writer(
                    fstream,
                    "dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r);" % Corr_s,
                )

    cw.writer(fstream, "// update wdot")
    for k in sorted(all_dict.keys()):
        s, nu = all_dict[k]
        if nu == 1:
            cw.writer(fstream, "wdot[%d] += q; // %s" % (k, s))
        elif nu == -1:
            cw.writer(fstream, "wdot[%d] -= q; // %s" % (k, s))
        elif nu > 0:
            cw.writer(fstream, "wdot[%d] += %.15g * q; // %s" % (k, nu, s))
        elif nu < 0:
            cw.writer(fstream, "wdot[%d] -= %.15g * q; // %s" % (k, -nu, s))

    if isPD:
        cw.writer(fstream, "// for convenience")
        cw.writer(fstream, "k_f *= Corr;")
        if reaction.reversible:
            cw.writer(fstream, "k_r *= Corr;")
    elif has_alpha:
        cw.writer(fstream, "// for convenience")
        cw.writer(fstream, "k_f *= alpha;")
        if reaction.reversible:
            cw.writer(fstream, "k_r *= alpha;")
        else:
            cw.writer(fstream, "k_r = 0.0;")

    if isPD:
        cw.writer(fstream, "dcdc_fac = 0.0;")
    # elif has_alpha:
    #    cw.writer(fstream,'dcdc_fac = q_nocor;')

    if has_alpha or isPD:

        # cw.writer(fstream,'if (consP) {')

        # for k in range(nSpecies):
        #    dqdc_s = self._Denhancement(mechanism,reaction,k,True)
        #    if dqdc_s != "0":
        #        if isPD:
        #            if dqdc_s == "1":
        #                dqdc_s ='dcdc_fac'
        #            else:
        #                dqdc_s +='*dcdc_fac'
        #        elif has_alpha:
        #            if dqdc_s == "1":
        #                dqdc_s ='q_nocor'
        #            else:
        #                dqdc_s +='*q_nocor'

        #    dqdc_s = dqdc_simple_precond(mechanism, species_info,reaction,sorted_reactants, sorted_products,  rea_dict, pro_dict, dqdc_s,k, removeForward)
        #    if dqdc_s:
        #        symb_k = species_info.nonqss_species[k].symbol
        #        cw.writer(fstream,'// d()/d[%s]' % symb_k)
        #        cw.writer(fstream,'dqdci = %s;' % (dqdc_s))
        #        #
        #        for m in sorted(all_dict.keys()):
        #            if all_dict[m][1] != 0:
        #                s1 = 'J[%d] += %.15g * dqdci;' % (k*(nSpecies+1)+m, all_dict[m][1])
        #                s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
        #                s2 = '// dwdot[%s]/d[%s]' % (all_dict[m][0], symb_k)
        #                cw.writer(fstream,s1.ljust(30) + s2)

        # cw.writer(fstream,'}')
        # cw.writer(fstream,'else {')

        for k in range(nSpecies):
            dqdc_s = Denhancement_d(
                mechanism, species_info, reaction, k, False
            )
            if dqdc_s != "0":
                if isPD:
                    if dqdc_s == "1":
                        dqdc_s = "dcdc_fac"
                    else:
                        dqdc_s += "*dcdc_fac"
                elif has_alpha:
                    if dqdc_s == "1":
                        dqdc_s = "q_nocor"
                    else:
                        dqdc_s += "*q_nocor"

            dqdc_s = dqdc_simple_precond(
                mechanism,
                species_info,
                reaction,
                sorted_reactants,
                sorted_products,
                rea_dict,
                pro_dict,
                dqdc_s,
                k,
                removeForward,
            )
            if dqdc_s:
                cw.writer(fstream, "dqdc[%d] = %s;" % (k, dqdc_s))
            else:
                cw.writer(fstream, "dqdc[%d] = 0.0;" % k)

        cw.writer(fstream, "for (int k=0; k<%d; k++) {" % nSpecies)
        for m in sorted(all_dict.keys()):
            if all_dict[m][1] != 0:
                s1 = "J[%d*k+%d] += %.15g * dqdc[k];" % (
                    (nSpecies + 1),
                    m,
                    all_dict[m][1],
                )
                s1 = s1.replace("+= 1 *", "+=").replace("+= -1 *", "-=")
                cw.writer(fstream, s1)
        # cw.writer(fstream,'}')

        cw.writer(fstream, "}")

        for m in sorted(all_dict.keys()):
            if all_dict[m][1] != 0:
                s1 = "J[%d] += %.15g * dqdT; // dwdot[%s]/dT" % (
                    nSpecies * (nSpecies + 1) + m,
                    all_dict[m][1],
                    all_dict[m][0],
                )
                s1 = s1.replace("+= 1 *", "+=").replace("+= -1 *", "-=")
                cw.writer(fstream, s1)

    else:

        for k in range(nSpecies):
            dqdc_s = dqdc_simple_precond(
                mechanism,
                species_info,
                reaction,
                sorted_reactants,
                sorted_products,
                rea_dict,
                pro_dict,
                "",
                k,
                removeForward,
            )
            if dqdc_s:
                cw.writer(fstream, "// d()/d[%s]" % all_dict[k][0])
                cw.writer(fstream, "dqdci = %s;" % (dqdc_s))
                if reaction.reversible or k in rea_dict:
                    for m in sorted(all_dict.keys()):
                        if all_dict[m][1] != 0:
                            s1 = "J[%d] += %.15g * dqdci;" % (
                                k * (nSpecies + 1) + m,
                                all_dict[m][1],
                            )
                            s1 = s1.replace("+= 1 *", "+=").replace(
                                "+= -1 *", "-="
                            )
                            s2 = "// dwdot[%s]/d[%s]" % (
                                all_dict[m][0],
                                all_dict[k][0],
                            )
                            cw.writer(fstream, s1.ljust(30) + s2)
        cw.writer(fstream, "// d()/dT")
        for m in sorted(all_dict.keys()):
            if all_dict[m][1] != 0:
                s1 = "J[%d] += %.15g * dqdT;" % (
                    nSpecies * (nSpecies + 1) + m,
                    all_dict[m][1],
                )
                s1 = (
                    s1.replace("+= 1 *", "+=")
                    .replace("+= -1 *", "-=")
                    .replace("+= -1 *", "-=")
                )
                s2 = "// dwdot[%s]/dT" % (all_dict[m][0])
                cw.writer(fstream, s1.ljust(30) + s2)


def ajac(fstream, mechanism, species_info, reaction_info):
    nSpecies = species_info.nSpecies
    nReactions = reaction_info.nReactions

    cw.writer(
        fstream,
    )
    cw.writer(fstream, cw.comment("compute the reaction Jacobian on GPU"))
    cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
    cw.writer(
        fstream,
        "void aJacobian(amrex::Real * J, amrex::Real * sc, amrex::Real T,"
        " const int consP)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
    )
    # Analytical jacobian not ready with QSS
    if species_info.nQSSspecies > 0:
        cw.writer(fstream, "// Do not use Analytical Jacobian with QSSA")
        cw.writer(fstream, "amrex::Abort();")
        cw.writer(
            fstream,
        )

    cw.writer(fstream, "for (int i=0; i<%d; i++) {" % (nSpecies + 1) ** 2)
    cw.writer(fstream, "J[i] = 0.0;")
    cw.writer(fstream, "}")

    cw.writer(
        fstream,
    )

    cw.writer(fstream, "amrex::Real wdot[%d];" % (nSpecies))
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies))
    cw.writer(fstream, "wdot[k] = 0.0;")
    cw.writer(fstream, "}")

    cw.writer(
        fstream,
    )

    cw.writer(
        fstream,
        "const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T }; //"
        " temperature cache",
    )
    cw.writer(fstream, "amrex::Real invT = 1.0 / tc[1];")
    cw.writer(fstream, "amrex::Real invT2 = invT * invT;")

    cw.writer(
        fstream,
    )

    if species_info.nQSSspecies > 0:
        cw.writer(fstream, "// Fill sc_qss here")
        cw.writer(
            fstream, "amrex::Real sc_qss[%d];" % species_info.nQSSspecies
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
        cw.writer(fstream, "comp_k_f_qss(tc, invT, kf_qss);")
        cw.writer(fstream, "comp_sc_qss(sc_qss, qf_qss, qr_qss);")
        cw.writer(
            fstream,
        )

    cw.writer(
        fstream,
    )

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

    cw.writer(
        fstream,
    )

    cw.writer(fstream, cw.comment("compute the mixture concentration"))
    cw.writer(fstream, "amrex::Real mixture = 0.0;")
    cw.writer(fstream, "for (int k = 0; k < %d; ++k) {" % nSpecies)
    cw.writer(fstream, "mixture += sc[k];")
    cw.writer(fstream, "}")

    cw.writer(
        fstream,
    )

    cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
    cw.writer(fstream, "amrex::Real g_RT[%d];" % (nSpecies))
    cw.writer(fstream, "gibbs(g_RT, tc);")
    if species_info.nQSSspecies > 0:
        cw.writer(
            fstream, "amrex::Real g_RT_qss[%d];" % (species_info.nQSSspecies)
        )
        cw.writer(fstream, "gibbs_qss(g_RT_qss, tc);")

    cw.writer(
        fstream,
    )

    cw.writer(fstream, cw.comment("compute the species enthalpy"))
    cw.writer(fstream, "amrex::Real h_RT[%d];" % (nSpecies))
    cw.writer(fstream, "speciesEnthalpy(h_RT, tc);")
    if species_info.nQSSspecies > 0:
        cw.writer(
            fstream, "amrex::Real h_RT_qss[%d];" % (species_info.nQSSspecies)
        )
        cw.writer(fstream, "speciesEnthalpy_qss(h_RT_qss, tc);")

    if species_info.nQSSspecies > 0:
        cw.writer(
            fstream,
        )
        cw.writer(fstream, "// Fill qss coeff")
        cw.writer(
            fstream,
            "comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, tc, g_RT, g_RT_qss);",
        )

    cw.writer(
        fstream,
    )

    cw.writer(
        fstream,
        "amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;",
    )
    cw.writer(fstream, "amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;")
    cw.writer(fstream, "amrex::Real dqdci, dcdc_fac, dqdc[%d];" % (nSpecies))
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

    for orig_idx, idx in reaction_info.idxmap.items():
        reaction = mechanism.reaction(orig_idx)
        # lt = reaction.lt
        # if lt:
        #     print("Landau-Teller reactions are not supported")
        #     sys.exit(1)

        cw.writer(
            fstream,
            cw.comment("reaction %d: %s" % (orig_idx, reaction.equation)),
        )
        falloff = reaction.reaction_type == "falloff"
        thirdBody = reaction.reaction_type == "three-body"
        if falloff:  # case 1
            cw.writer(fstream, cw.comment("a pressure-fall-off reaction"))
            ajac_reaction_d(
                fstream,
                mechanism,
                species_info,
                reaction_info,
                reaction,
                orig_idx,
                1,
            )
        elif thirdBody:  # case 2
            cw.writer(
                fstream,
                cw.comment("a third-body and non-pressure-fall-off reaction"),
            )
            ajac_reaction_d(
                fstream,
                mechanism,
                species_info,
                reaction_info,
                reaction,
                orig_idx,
                2,
            )
        else:  # case 3
            cw.writer(
                fstream,
                cw.comment(
                    "a non-third-body and non-pressure-fall-off reaction"
                ),
            )
            ajac_reaction_d(
                fstream,
                mechanism,
                species_info,
                reaction_info,
                reaction,
                orig_idx,
                3,
            )
        cw.writer(
            fstream,
        )

    cw.writer(
        fstream,
        "amrex::Real c_R[%d], dcRdT[%d], e_RT[%d];"
        % (nSpecies, nSpecies, nSpecies),
    )
    cw.writer(fstream, "amrex::Real * eh_RT;")
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

    cw.writer(
        fstream,
    )

    cw.writer(
        fstream,
        "amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;",
    )
    cw.writer(fstream, "for (int k = 0; k < %d; ++k) {" % nSpecies)
    cw.writer(fstream, "cmix += c_R[k]*sc[k];")
    cw.writer(fstream, "dcmixdT += dcRdT[k]*sc[k];")
    cw.writer(fstream, "ehmix += eh_RT[k]*wdot[k];")
    cw.writer(
        fstream,
        "dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[%d+k];"
        % (nSpecies * (nSpecies + 1)),
    )
    cw.writer(fstream, "}")

    cw.writer(
        fstream,
    )
    cw.writer(fstream, "amrex::Real cmixinv = 1.0/cmix;")
    cw.writer(fstream, "amrex::Real tmp1 = ehmix*cmixinv;")
    cw.writer(fstream, "amrex::Real tmp3 = cmixinv*T;")
    cw.writer(fstream, "amrex::Real tmp2 = tmp1*tmp3;")
    cw.writer(fstream, "amrex::Real dehmixdc;")

    cw.writer(fstream, "// dTdot/d[X]")
    cw.writer(fstream, "for (int k = 0; k < %d; ++k) {" % nSpecies)
    cw.writer(fstream, "dehmixdc = 0.0;")
    cw.writer(fstream, "for (int m = 0; m < %d; ++m) {" % nSpecies)
    cw.writer(fstream, "dehmixdc += eh_RT[m]*J[k*%s+m];" % (nSpecies + 1))
    cw.writer(fstream, "}")
    cw.writer(
        fstream,
        "J[k*%d+%d] = tmp2*c_R[k] - tmp3*dehmixdc;" % (nSpecies + 1, nSpecies),
    )
    cw.writer(fstream, "}")

    cw.writer(fstream, "// dTdot/dT")
    cw.writer(
        fstream,
        "J[%d] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;"
        % (nSpecies * (nSpecies + 1) + nSpecies),
    )

    cw.writer(
        fstream,
    )
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def ajac_reaction_d(
    fstream, mechanism, species_info, reaction_info, reaction, orig_idx, rcase
):

    nSpecies = species_info.nSpecies
    removeForward = cu.isRemoveForward(reaction_info, orig_idx)
    thirdBody = reaction.reaction_type == "three-body"
    if rcase == 1:  # pressure-dependent reaction
        isPD = True
        # FIXME is this right?
        has_alpha = True
        cw.writer(fstream, "// also 3-body")
        # if reaction.thirdBody:
        #     has_alpha = True
        #     cw.writer(fstream,"// also 3-body")
        # else:
        #     has_alpha = False
        #     cw.writer(fstream,"// non 3-body")
        #     print(
        #         "FIXME: pressure dependent non-3-body reaction in _ajac_reaction"
        #     )
        #     sys.exit(1)
    elif rcase == 2:  # third-body and non-pressure-dependent reaction
        isPD = False
        has_alpha = True
    elif rcase == 3:  # simple non-third and non-pressure-dependent reaction
        isPD = False
        has_alpha = False
    else:
        print("_ajac_reaction: wrong case ", rcase)
        exit(1)

    rea_dict = OrderedDict()
    pro_dict = OrderedDict()
    all_dict = OrderedDict()
    sumNuk = 0
    for symbol, coefficient in reaction.reactants.items():
        k = species_info.ordered_idx_map[symbol]
        sumNuk -= coefficient
        if k in rea_dict:
            coe_old = rea_dict[k][1]
            rea_dict[k] = (symbol, coefficient + coe_old)
        else:
            rea_dict[k] = (symbol, coefficient)
    for symbol, coefficient in reaction.products.items():
        k = species_info.ordered_idx_map[symbol]
        sumNuk += coefficient
        if k in pro_dict:
            coe_old = pro_dict[k][1]
            pro_dict[k] = (symbol, coefficient + coe_old)
        else:
            pro_dict[k] = (symbol, coefficient)
    for k in range(nSpecies):
        # QSS at the end so we should be good
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

    sorted_reactants = sorted(rea_dict.values())
    sorted_products = sorted(pro_dict.values())

    if not reaction.reversible:
        if isPD or has_alpha:
            print(
                "FIXME: irreversible reaction in _ajac_reaction may not work"
            )
            cw.writer(
                fstream,
                "// FIXME: irreversible reaction in _ajac_reaction may not"
                " work",
            )
        for k in range(nSpecies):
            if k in sorted_reactants and k in sorted_products:
                print(
                    "FIXME: irreversible reaction in _ajac_reaction may not"
                    " work"
                )
                cw.writer(
                    fstream,
                    "// FIXME: irreversible reaction in _ajac_reaction may not"
                    " work",
                )

    dim = cu.phaseSpaceUnits(reaction.reactants)
    aeuc = cu.activationEnergyUnits()
    is_sri = False
    is_troe = False
    if isPD:
        Corr_s = "Corr *"
        uc = cu.prefactorUnits(
            cc.ureg("mole/cm**3"), 1 - dim
        )  # Case 2 !PD, TB
        ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), 1 - dim)
        A = (reaction.high_rate.pre_exponential_factor * ctuc).to(
            uc.to_root_units()
        )
        beta = reaction.high_rate.temperature_exponent
        E = (
            reaction.high_rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
        ).to(aeuc)

        low_A = (reaction.low_rate.pre_exponential_factor * ctuc).to(
            uc.to_root_units()
        )
        low_beta = reaction.low_rate.temperature_exponent
        low_E = (
            reaction.low_rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
        ).to(aeuc)
        if reaction.rate.type == "Troe":
            troe = reaction.rate.falloff_coeffs
            ntroe = len(troe)
            is_troe = True
        elif reaction.rate.type == "Sri":
            sri = reaction.rate.falloff_coeffs
            nsri = len(sri)
            is_sri = True
        else:
            print("Unrecognized reaction rate type")
            sys.exit(1)
    elif has_alpha:
        Corr_s = "alpha * "
        uc = cu.prefactorUnits(cc.ureg("mole/cm**3"), -dim)  # Case 2 !PD, TB
        ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), -dim)
        A = (reaction.rate.pre_exponential_factor * ctuc).to(
            uc.to_root_units()
        )
        beta = reaction.rate.temperature_exponent
        E = (
            reaction.rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
        ).to(aeuc)
    else:
        uc = cu.prefactorUnits(
            cc.ureg("mole/cm**3"), 1 - dim
        )  # Case 3 !PD, !TB
        ctuc = cu.prefactorUnits(cc.ureg("kmol/m**3"), 1 - dim)
        A = (reaction.rate.pre_exponential_factor * ctuc).to(
            uc.to_root_units()
        )
        beta = reaction.rate.temperature_exponent
        E = (
            reaction.rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
        ).to(aeuc)
        Corr_s = ""

    if has_alpha:
        cw.writer(fstream, "// 3-body correction factor")
        cw.writer(
            fstream,
            "alpha = %s;" % enhancement_d(mechanism, species_info, reaction),
        )

    # forward
    cw.writer(fstream, "// forward")
    cw.writer(
        fstream,
        "phi_f = %s;"
        % cu.QSSsortedPhaseSpace(mechanism, species_info, reaction.reactants),
    )
    #
    cw.writer(fstream, "k_f = %.15g" % (A.m))
    cw.writer(
        fstream,
        "            * exp(%.15g * tc[0] - %.15g * (%.15g) * invT);"
        % (beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, E.m),
    )
    #
    if removeForward:
        cw.writer(fstream, "// Remove forward reaction")
        # DLNKFDT CHECK
        cw.writer(
            fstream,
            "dlnkfdT = %.15g * invT + %.15g *  %.15g  * invT2;"
            % (beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, E.m),
        )
        cw.writer(fstream, "//dlnkfdT = 0.0;")
    else:
        cw.writer(
            fstream,
            "dlnkfdT = %.15g * invT + %.15g *  %.15g  * invT2;"
            % (beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, E.m),
        )

    if isPD:
        cw.writer(fstream, "// pressure-fall-off")
        cw.writer(
            fstream,
            "k_0 = %.15g * exp(%.15g * tc[0] - %.15g * (%.15g) * invT);"
            % (
                low_A.m * 10 ** (3**dim),
                low_beta,
                (1.0 / cc.Rc / cc.ureg.kelvin).m,
                low_E.m,
            ),
        )
        cw.writer(fstream, "Pr = 1e-%d * alpha / k_f * k_0;" % (dim * 6))
        cw.writer(fstream, "fPr = Pr / (1.0+Pr);")
        cw.writer(
            fstream,
            "dlnk0dT = %.15g * invT + %.15g * (%.15g) * invT2;"
            % (low_beta, (1.0 / cc.Rc / cc.ureg.kelvin).m, low_E.m),
        )
        cw.writer(fstream, "dlogPrdT = log10e*(dlnk0dT - dlnkfdT);")
        cw.writer(fstream, "dlogfPrdT = dlogPrdT / (1.0+Pr);")
        #
        if is_sri:
            cw.writer(fstream, "// SRI form")
            print("FIXME: sri not supported in _ajac_reaction yet")
            sys.exit(1)
        elif is_troe:
            cw.writer(fstream, "// Troe form")
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
            cw.writer(fstream, "// Lindemann form")
            cw.writer(fstream, "F = 1.0;")
            cw.writer(fstream, "dlogFdlogPr = 0.0;")
            cw.writer(fstream, "dlogFdT = 0.0;")

    # reverse
    if not reaction.reversible:
        cw.writer(fstream, "// rate of progress")
        if (not has_alpha) and (not isPD):
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(fstream, "//q = k_f*phi_f;")
                cw.writer(fstream, "q = 0;")
            else:
                cw.writer(fstream, "q = k_f*phi_f;")
        else:
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(fstream, "//q_nocor = k_f*phi_f;")
                cw.writer(fstream, "q_nocor = 0;")
            else:
                cw.writer(fstream, "q_nocor = k_f*phi_f;")

            if isPD:
                cw.writer(fstream, "Corr = fPr * F;")
                cw.writer(fstream, "q = Corr * q_nocor;")
            else:
                cw.writer(fstream, "q = alpha * q_nocor;")

        if isPD:
            cw.writer(fstream, "dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);")
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(
                    fstream,
                    "//dqdT = %sdlnkfdT*k_f*phi_f + dlnCorrdT*q;" % Corr_s,
                )
                cw.writer(fstream, "dqdT =  dlnCorrdT*q;")
            else:
                cw.writer(
                    fstream,
                    "dqdT = %sdlnkfdT*k_f*phi_f + dlnCorrdT*q;" % Corr_s,
                )
        else:
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(fstream, "//dqdT = %sdlnkfdT*k_f*phi_f;" % Corr_s)
                cw.writer(fstream, "dqdT = 0;")
            else:
                cw.writer(fstream, "dqdT = %sdlnkfdT*k_f*phi_f;" % Corr_s)
    else:
        cw.writer(fstream, "// reverse")
        cw.writer(
            fstream,
            "phi_r = %s;"
            % cu.QSSsortedPhaseSpace(
                mechanism, species_info, reaction.products
            ),
        )
        cw.writer(
            fstream,
            "Kc = %s;" % cu.sortedKc(mechanism, species_info, reaction),
        )
        cw.writer(fstream, "k_r = k_f / Kc;")

        dlnKcdT_s = "invT * ("
        terms = []
        for symbol, coefficient in sorted(
            sorted_reactants,
            key=lambda x: next(
                (y for y in species_info.all_species if y.name == x[0]),
                None,
            ).idx,
        ):
            if symbol not in species_info.qss_species_list:
                k = species_info.ordered_idx_map[symbol]
                if coefficient == 1.0:
                    terms.append("h_RT[%d]" % (k))
                else:
                    terms.append("%f*h_RT[%d]" % (coefficient, k))
            else:
                k = species_info.ordered_idx_map[symbol] - nSpecies
                if coefficient == 1.0:
                    terms.append("h_RT_qss[%d]" % (k))
                else:
                    terms.append("%f*h_RT_qss[%d]" % (coefficient, k))
        dlnKcdT_s += "-(" + " + ".join(terms) + ")"
        terms = []
        for symbol, coefficient in sorted(
            sorted_products,
            key=lambda x: next(
                (y for y in species_info.all_species if y.name == x[0]),
                None,
            ).idx,
        ):
            if symbol not in species_info.qss_species_list:
                k = species_info.ordered_idx_map[symbol]
                if coefficient == 1.0:
                    terms.append("h_RT[%d]" % (k))
                else:
                    terms.append("%f*h_RT[%d]" % (coefficient, k))
            else:
                k = species_info.ordered_idx_map[symbol] - nSpecies
                if coefficient == 1.0:
                    terms.append("h_RT_qss[%d]" % (k))
                else:
                    terms.append("%f*h_RT_qss[%d]" % (coefficient, k))
        dlnKcdT_s += " + (" + " + ".join(terms) + ")"
        if sumNuk > 0:
            dlnKcdT_s += " - %f" % sumNuk
        elif sumNuk < 0:
            dlnKcdT_s += " + %f" % (-sumNuk)
        dlnKcdT_s += ")"
        cw.writer(fstream, "dlnKcdT = %s;" % dlnKcdT_s)

        cw.writer(fstream, "dkrdT = (dlnkfdT - dlnKcdT)*k_r;")

        cw.writer(fstream, "// rate of progress")
        if (not has_alpha) and (not isPD):
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(fstream, "//q = k_f*phi_f - k_r*phi_r;")
                cw.writer(fstream, "q = - k_r*phi_r;")
            else:
                cw.writer(fstream, "q = k_f*phi_f - k_r*phi_r;")
        else:
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(fstream, "//q_nocor = k_f*phi_f - k_r*phi_r;")
                cw.writer(fstream, "q_nocor = - k_r*phi_r;")
            else:
                cw.writer(fstream, "q_nocor = k_f*phi_f - k_r*phi_r;")
            if isPD:
                cw.writer(fstream, "Corr = fPr * F;")
                cw.writer(fstream, "q = Corr * q_nocor;")
            else:
                cw.writer(fstream, "q = alpha * q_nocor;")

        if isPD:
            cw.writer(fstream, "dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);")
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(
                    fstream,
                    "//dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r) +"
                    " dlnCorrdT*q;" % Corr_s,
                )
                cw.writer(
                    fstream, "dqdT = %s(- dkrdT*phi_r) + dlnCorrdT*q;" % Corr_s
                )
            else:
                cw.writer(
                    fstream,
                    "dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;"
                    % Corr_s,
                )
        else:
            if removeForward:
                cw.writer(fstream, "// Remove forward reaction")
                cw.writer(
                    fstream,
                    "//dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r);" % Corr_s,
                )
                cw.writer(fstream, "dqdT = %s( - dkrdT*phi_r);" % Corr_s)
            else:
                cw.writer(
                    fstream,
                    "dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r);" % Corr_s,
                )

    cw.writer(fstream, "// update wdot")
    # only the nSpecies transported in all_dict
    for k in sorted(all_dict.keys()):
        s, nu = all_dict[k]
        if nu == 1:
            cw.writer(fstream, "wdot[%d] += q; // %s" % (k, s))
        elif nu == -1:
            cw.writer(fstream, "wdot[%d] -= q; // %s" % (k, s))
        elif nu > 0:
            cw.writer(fstream, "wdot[%d] += %.15g * q; // %s" % (k, nu, s))
        elif nu < 0:
            cw.writer(fstream, "wdot[%d] -= %.15g * q; // %s" % (k, -nu, s))

    if isPD:
        cw.writer(fstream, "// for convenience")
        cw.writer(fstream, "k_f *= Corr;")
        if reaction.reversible:
            cw.writer(fstream, "k_r *= Corr;")
    elif has_alpha:
        cw.writer(fstream, "// for convenience")
        cw.writer(fstream, "k_f *= alpha;")
        if reaction.reversible:
            cw.writer(fstream, "k_r *= alpha;")
        else:
            cw.writer(fstream, "k_r = 0.0;")

    if isPD:
        cw.writer(fstream, "dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);")
    # elif has_alpha:
    #    cw.writer(fstream,'dcdc_fac = q_nocor;')

    if has_alpha or isPD:

        cw.writer(fstream, "if (consP) {")

        for k in range(nSpecies):
            dqdc_s = Denhancement_d(mechanism, species_info, reaction, k, True)
            if dqdc_s != "0":
                if isPD:
                    if dqdc_s == "1":
                        dqdc_s = "dcdc_fac"
                    else:
                        dqdc_s += "*dcdc_fac"
                elif has_alpha:
                    if dqdc_s == "1":
                        dqdc_s = "q_nocor"
                    else:
                        dqdc_s += "*q_nocor"

            dqdc_s = dqdc_simple_d(
                mechanism,
                species_info,
                reaction,
                sorted_reactants,
                sorted_products,
                rea_dict,
                pro_dict,
                dqdc_s,
                k,
                removeForward,
            )
            if dqdc_s:
                symb_k = species_info.nonqss_species[k].name
                cw.writer(fstream, "// d()/d[%s]" % symb_k)
                cw.writer(fstream, "dqdci = %s;" % (dqdc_s))
                #
                for m in sorted(all_dict.keys()):
                    if all_dict[m][1] != 0:
                        s1 = "J[%d] += %.15g * dqdci;" % (
                            k * (nSpecies + 1) + m,
                            all_dict[m][1],
                        )
                        s1 = s1.replace("+= 1 *", "+=").replace(
                            "+= -1 *", "-="
                        )
                        s2 = "// dwdot[%s]/d[%s]" % (
                            all_dict[m][0],
                            symb_k,
                        )
                        cw.writer(fstream, s1.ljust(30) + s2)

        cw.writer(fstream, "}")
        cw.writer(fstream, "else {")

        for k in range(nSpecies):
            dqdc_s = Denhancement_d(
                mechanism, species_info, reaction, k, False
            )
            if dqdc_s != "0":
                if isPD:
                    if dqdc_s == "1":
                        dqdc_s = "dcdc_fac"
                    else:
                        dqdc_s += "*dcdc_fac"
                elif has_alpha:
                    if dqdc_s == "1":
                        dqdc_s = "q_nocor"
                    else:
                        dqdc_s += "*q_nocor"

            dqdc_s = dqdc_simple_d(
                mechanism,
                species_info,
                reaction,
                sorted_reactants,
                sorted_products,
                rea_dict,
                pro_dict,
                dqdc_s,
                k,
                removeForward,
            )
            if dqdc_s:
                cw.writer(fstream, "dqdc[%d] = %s;" % (k, dqdc_s))

        cw.writer(fstream, "for (int k=0; k<%d; k++) {" % nSpecies)
        for m in sorted(all_dict.keys()):
            if all_dict[m][1] != 0:
                s1 = "J[%d*k+%d] += %.15g * dqdc[k];" % (
                    (nSpecies + 1),
                    m,
                    all_dict[m][1],
                )
                s1 = s1.replace("+= 1 *", "+=").replace("+= -1 *", "-=")
                cw.writer(fstream, s1)
        cw.writer(fstream, "}")

        cw.writer(fstream, "}")

        for m in sorted(all_dict.keys()):
            if all_dict[m][1] != 0:
                s1 = "J[%d] += %.15g * dqdT; // dwdot[%s]/dT" % (
                    nSpecies * (nSpecies + 1) + m,
                    all_dict[m][1],
                    all_dict[m][0],
                )
                s1 = s1.replace("+= 1 *", "+=").replace("+= -1 *", "-=")
                cw.writer(fstream, s1)

    else:

        for k in range(nSpecies):
            dqdc_s = dqdc_simple_d(
                mechanism,
                species_info,
                reaction,
                sorted_reactants,
                sorted_products,
                rea_dict,
                pro_dict,
                "",
                k,
                removeForward,
            )
            if dqdc_s:
                cw.writer(fstream, "// d()/d[%s]" % all_dict[k][0])
                cw.writer(fstream, "dqdci = %s;" % (dqdc_s))
                if reaction.reversible or k in rea_dict:
                    for m in sorted(all_dict.keys()):
                        if all_dict[m][1] != 0:
                            s1 = "J[%d] += %.15g * dqdci;" % (
                                k * (nSpecies + 1) + m,
                                all_dict[m][1],
                            )
                            s1 = s1.replace("+= 1 *", "+=").replace(
                                "+= -1 *", "-="
                            )
                            s2 = "// dwdot[%s]/d[%s]" % (
                                all_dict[m][0],
                                all_dict[k][0],
                            )
                            cw.writer(fstream, s1.ljust(30) + s2)
        cw.writer(fstream, "// d()/dT")
        for m in sorted(all_dict.keys()):
            if all_dict[m][1] != 0:
                s1 = "J[%d] += %.15g * dqdT;" % (
                    nSpecies * (nSpecies + 1) + m,
                    all_dict[m][1],
                )
                s1 = (
                    s1.replace("+= 1 *", "+=")
                    .replace("+= -1 *", "-=")
                    .replace("+= -1 *", "-=")
                )
                s2 = "// dwdot[%s]/dT" % (all_dict[m][0])
                cw.writer(fstream, s1.ljust(30) + s2)


def dqdc_simple_precond(
    mechanism,
    species_info,
    reaction,
    sorted_reactants,
    sorted_products,
    rea_dict,
    pro_dict,
    dqdc_s,
    k,
    removeForward,
):
    if dqdc_s == "0":
        dqdc_s = ""
    if k in sorted(rea_dict.keys()):
        dps = DphaseSpace(
            mechanism, species_info, sorted_reactants, rea_dict[k][0]
        )
        if dps == "1.0":
            dps_s = ""
        else:
            dps_s = "*" + dps
        if removeForward:
            cw.writer(fstream, "// Remove forward reaction")
            dqdc_s += ""
        else:
            dqdc_s += " + k_f%s" % dps_s
    if reaction.reversible:
        if k in sorted(pro_dict.keys()):
            dps = DphaseSpace(
                mechanism, species_info, sorted_products, pro_dict[k][0]
            )
            if dps == "1.0":
                dps_s = ""
            else:
                dps_s = "*" + dps
            dqdc_s += " - k_r%s" % dps_s
    return dqdc_s


def dqdc_simple_d(
    mechanism,
    species_info,
    reaction,
    sorted_reactants,
    sorted_products,
    rea_dict,
    pro_dict,
    dqdc_s,
    k,
    removeForward,
):
    if dqdc_s == "0":
        dqdc_s = ""
    if k in sorted(rea_dict.keys()):
        dps = DphaseSpace(
            mechanism, species_info, sorted_reactants, rea_dict[k][0]
        )
        if dps == "1.0":
            dps_s = ""
        else:
            dps_s = "*" + dps
        if removeForward:
            cw.writer(fstream, "// Remove forward reaction")
            dqdc_s += ""
        else:
            dqdc_s += " + k_f%s" % dps_s
    if reaction.reversible:
        if k in sorted(pro_dict.keys()):
            dps = DphaseSpace(
                mechanism, species_info, sorted_products, pro_dict[k][0]
            )
            if dps == "1.0":
                dps_s = ""
            else:
                dps_s = "*" + dps
            dqdc_s += " - k_r%s" % dps_s
    return dqdc_s


def enhancement_d(mechanism, species_info, reaction):
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
    for i, (symbol, efficiency) in enumerate(efficiencies.items()):
        if symbol not in species_info.qss_species_list:
            factor = "( %.15g - 1)" % (efficiency)
            conc = "sc[%d]" % species_info.ordered_idx_map[symbol]
            alpha.append("%s*%s" % (factor, conc))

    return " + ".join(alpha).replace("+ -", "- ")


def Denhancement_d(mechanism, species_info, reaction, kid, consP):
    thirdBody = reaction.reaction_type == "three-body"
    falloff = reaction.reaction_type == "falloff"
    if not thirdBody and not falloff:
        print("Denhancement_d called for a reaction without a third body")
        sys.exit(1)

    if not hasattr(reaction, "efficiencies"):
        print("FIXME")
        sys.exit(1)
        species, coefficient = thirdBody
        if species == "<mixture>":
            if consP:
                return "0"
            else:
                return "1"
        elif species_info.ordered_idx_map[species] == kid:
            return "1"
        else:
            return "0"
    else:
        efficiencies = reaction.efficiencies
        if consP:
            for i, (symbol, efficiency) in enumerate(efficiencies.items()):
                if species_info.ordered_idx_map[symbol] == kid:
                    return "(%.15g - 1)" % (efficiency)
            return "0"
        else:
            for i, (symbol, efficiency) in enumerate(efficiencies.items()):
                if species_info.ordered_idx_map[symbol] == kid:
                    return "%.15g" % (efficiency)
            return "1"


def DphaseSpace(mechanism, species_info, reagents, r):
    phi = []

    for symbol, coefficient in sorted(
        reagents,
        key=lambda x: next(
            (y for y in species_info.all_species if y.name == x[0]),
            None,
        ).idx,
    ):
        if symbol not in species_info.qss_species_list:
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
        else:
            if symbol == r:
                if coefficient > 1:
                    phi += ["%f" % coefficient]
                    if (coefficient - 1) == 1.0:
                        conc = "sc_qss[%d]" % (
                            species_info.ordered_idx_map[symbol]
                            - species_info.nSpecies
                        )
                    else:
                        conc = "pow(sc_qss[%d],%f)" % (
                            species_info.ordered_idx_map[symbol]
                            - species_info.nSpecies,
                            (coefficient - 1),
                        )
                    phi += [conc]
            else:
                if coefficient == 1.0:
                    conc = "sc_qss[%d]" % (
                        species_info.ordered_idx_map[symbol]
                        - species_info.nSpecies
                    )
                else:
                    conc = "pow(sc_qss[%d], %f)" % (
                        species_info.ordered_idx_map[symbol]
                        - species_info.nSpecies,
                        coefficient,
                    )
                phi += [conc]

    if phi:
        return "*".join(phi)
    else:
        return "1.0"


def DproductionRatePrecond(fstream, mechanism, species_info, reaction_info):
    nSpecies = species_info.nSpecies
    cw.writer(fstream)
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
    cw.writer(fstream, "{")

    cw.writer(fstream, "amrex::Real c[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % nSpecies)
    cw.writer(fstream, "c[k] = 1.e6 * sc[k];")
    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "aJacobian_precond(J, c, *Tp, *HP);")

    cw.writer(fstream)
    cw.writer(fstream, "// dwdot[k]/dT")
    cw.writer(fstream, "// dTdot/d[X]")
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % nSpecies)
    cw.writer(fstream, "J[%d+k] *= 1.e-6;" % (nSpecies * (nSpecies + 1)))
    cw.writer(fstream, "J[k*%d+%d] *= 1.e6;" % (nSpecies + 1, nSpecies))
    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")


def DproductionRate(fstream, mechanism, species_info, reaction_info):
    nSpecies = species_info.nSpecies

    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the reaction Jacobian"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *  J,"
        " amrex::Real *  sc, amrex::Real *  Tp, const int * consP)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "amrex::Real c[%d];" % (nSpecies))
    cw.writer(fstream)
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % nSpecies)
    cw.writer(fstream, "c[k] = 1.e6 * sc[k];")
    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "aJacobian(J, c, *Tp, *consP);")

    cw.writer(fstream)
    cw.writer(fstream, "// dwdot[k]/dT")
    cw.writer(fstream, "// dTdot/d[X]")
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % nSpecies)
    cw.writer(fstream, "J[%d+k] *= 1.e-6;" % (nSpecies * (nSpecies + 1)))
    cw.writer(fstream, "J[k*%d+%d] *= 1.e6;" % (nSpecies + 1, nSpecies))
    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")
    cw.writer(fstream)