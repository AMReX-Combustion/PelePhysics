"""Symbolic math for symbolic differentiation."""
import re
import time
from collections import OrderedDict

import numpy as np
import symengine as sme
import sympy as smp

import ceptr.thermo as cth


class SymbolicMath:
    """Symbols to carry throughout operations."""

    def __init__(self, species_info, reaction_info, mechanism):

        n_species = species_info.n_species
        n_qssa_species = species_info.n_qssa_species

        self.T_smp = smp.symbols("T")
        self.tc_smp = [
            smp.log(self.T_smp),
            self.T_smp,
            self.T_smp**2,
            self.T_smp**3,
            self.T_smp**4,
        ]
        self.invT_smp = 1.0 / self.tc_smp[1]
        self.invT2_smp = self.invT_smp * self.invT_smp

        self.sc_smp = [
            smp.symbols("sc[" + str(i) + "]") for i in range(n_species)
        ]
        self.g_RT_smp = [
            smp.symbols("g_RT[" + str(i) + "]") for i in range(n_species)
        ]
        self.h_RT_smp = [
            smp.symbols("h_RT[" + str(i) + "]") for i in range(n_species)
        ]
        self.wdot_smp = [
            smp.symbols("wdot[" + str(i) + "]") for i in range(n_species)
        ]
        self.jac_smp = [
            smp.symbols(f"J[{str(i)}][{str(j)}]")
            for i in range(n_species)
            for j in range(n_species + 1)
        ]

        # temporary symbols that contain temperature dependence

        species_coeffs = cth.analyze_thermodynamics(mechanism, species_info, 0)
        low_temp, high_temp, midpoints = species_coeffs
        self.midpointsList = []
        for mid_temp, species_list in list(midpoints.items()):
            self.midpointsList.append(mid_temp)
        self.midpointsList_sorted = sorted(self.midpointsList)

        self.g_RT_smp_tmp = {}
        self.h_RT_smp_tmp = {}

        for mid_temp in self.midpointsList_sorted:
            self.g_RT_smp_tmp[mid_temp] = {}
            self.h_RT_smp_tmp[mid_temp] = {}
            self.g_RT_smp_tmp[mid_temp]["m"] = [
                smp.symbols("g_RT[" + str(i) + "]") for i in range(n_species)
            ]
            self.g_RT_smp_tmp[mid_temp]["p"] = [
                smp.symbols("g_RT[" + str(i) + "]") for i in range(n_species)
            ]
            self.h_RT_smp_tmp[mid_temp]["m"] = [
                smp.symbols("h_RT[" + str(i) + "]") for i in range(n_species)
            ]
            self.h_RT_smp_tmp[mid_temp]["p"] = [
                smp.symbols("h_RT[" + str(i) + "]") for i in range(n_species)
            ]

        # mixture (useful for third body reactions)
        self.mixture_smp = 0.0
        for i in range(n_species):
            self.mixture_smp += self.sc_smp[i]

        # for storing intermediate terms used in definition of sc_qss
        self.intermediate_helpers_smp = OrderedDict()
        self.intermediate_terms_smp = {}

        if n_qssa_species > 0:

            n_qssa_reactions = reaction_info.n_qssa_reactions

            self.sc_qss_smp = [
                smp.symbols("sc_qss[" + str(i) + "]")
                for i in range(n_qssa_species)
            ]
            self.kf_qss_smp = [
                smp.symbols("kf_qss[" + str(i) + "]")
                for i in range(n_qssa_reactions)
            ]
            self.qf_qss_smp = [
                smp.symbols("qf_qss[" + str(i) + "]")
                for i in range(n_qssa_reactions)
            ]
            self.qr_qss_smp = [
                smp.symbols("qr_qss[" + str(i) + "]")
                for i in range(n_qssa_reactions)
            ]
            self.g_RT_qss_smp = [
                smp.symbols("g_RT_qss[" + str(i) + "]")
                for i in range(n_qssa_species)
            ]
            self.h_RT_qss_smp = [
                smp.symbols("h_RT_qss[" + str(i) + "]")
                for i in range(n_qssa_species)
            ]

            # temporary symbols that contain temperature dependence
            qss_species_coeffs = cth.analyze_thermodynamics(
                mechanism, species_info, 1
            )
            low_temp, high_temp, midpoints = qss_species_coeffs
            self.midpointsQSSList = []
            for mid_temp, species_list in list(midpoints.items()):
                self.midpointsQSSList.append(mid_temp)
            self.midpointsQSSList_sorted = sorted(self.midpointsQSSList)

            self.g_RT_qss_smp_tmp = {}
            self.h_RT_qss_smp_tmp = {}

            for mid_temp in self.midpointsQSSList_sorted:
                self.g_RT_qss_smp_tmp[mid_temp] = {}
                self.h_RT_qss_smp_tmp[mid_temp] = {}
                self.g_RT_qss_smp_tmp[mid_temp]["m"] = [
                    smp.symbols("g_RT_qss[" + str(i) + "]")
                    for i in range(n_qssa_species)
                ]
                self.g_RT_qss_smp_tmp[mid_temp]["p"] = [
                    smp.symbols("g_RT_qss[" + str(i) + "]")
                    for i in range(n_qssa_species)
                ]
                self.h_RT_qss_smp_tmp[mid_temp]["m"] = [
                    smp.symbols("h_RT_qss[" + str(i) + "]")
                    for i in range(n_qssa_species)
                ]
                self.h_RT_qss_smp_tmp[mid_temp]["p"] = [
                    smp.symbols("h_RT_qss[" + str(i) + "]")
                    for i in range(n_qssa_species)
                ]

            self.kf_qss_smp_tmp = [
                smp.symbols("kf_qss[" + str(i) + "]")
                for i in range(n_qssa_reactions)
            ]

            # Create dict to hold end of chain rule dscqssdsc terms
            self.dscqssdsc_stop = {"info": ""}

    def convert_to_cpp(self, sym_smp):
        """Convert sympy object to C code compatible string."""
        # Convert to ccode (to fix pow) and then string
        cppcode = smp.ccode(sym_smp)
        cpp_str = str(cppcode)

        return cpp_str

    def write_array_to_cpp(self, list_smp, array_str, cw, fstream):
        """Convert sympy array to C code compatible string."""
        n = len(list_smp)

        # Write common expressions
        times = time.time()
        array_cse = smp.cse(list_smp)
        for cse_idx in range(len(array_cse[0])):
            left_cse = self.convert_to_cpp(array_cse[0][cse_idx][0])
            right_cse = self.convert_to_cpp(array_cse[0][cse_idx][1])
            cw.writer(
                fstream,
                "const amrex::Real %s = %s;"
                % (
                    left_cse,
                    right_cse,
                ),
            )
        timee = time.time()

        print("Made common expr (time = %.3g s)" % (timee - times))

        # Write all the entries
        for i in range(n):
            # The full expression is stored in array_cse index 1
            times = time.time()
            cpp_str = self.convert_to_cpp(array_cse[1][i])
            timee = time.time()
            print(
                "Made expr for entry %d (time = %.3g s)" % (i, timee - times)
            )
            times = time.time()
            cw.writer(
                fstream,
                "%s[%s] = %s;"
                % (
                    array_str,
                    str(i),
                    cpp_str,
                ),
            )
            timee = time.time()
            print(
                "Printed expr for entry %d (time = %.3g s)"
                % (i, timee - times)
            )

    def syms_to_specnum(self, sym_smp):
        """Extracts number from syms string"""
        num = re.findall(r"\[(.*?)\]", str(sym_smp))
        return int(num[0])

    def chain_diff(self, species_info, sym_smp, ref_var):
        """Chain rule diff that replacing sc_qss with dependents."""

        free_symb = sym_smp.free_symbols
        # create dict of sc_terms as function of sc_qss terms
        sc_terms = {}
        for sc_symb in free_symb:
            if "sc_qss" in str(sc_symb):
                scqssnum = self.syms_to_specnum(sc_symb)
                for name in species_info.dict_qssdepend_sc[
                    species_info.qssa_species_list[scqssnum]
                ]:
                    if name in sc_terms:
                        sc_terms[name].append(sc_symb)
                    else:
                        sc_terms[name] = [sc_symb]

        # # Let's try subs
        # print(f"Substituting sc_qss terms...")
        # if ref_var in sc_terms:
        #     for scqss in sc_terms[ref_var]:
        #         print(scqss)
        #         scqssnum = self.syms_to_specnum(scqss)
        #         sym_smp = sym_smp.subs(self.sc_qss_smp[scqssnum], scqss)

        # print(sym_smp.free_symbols)
        # exit()

        print(f"Computing dsym/d{ref_var}...")
        dsym_dref = sme.diff(sym_smp, ref_var)

        # check for sc_qss terms
        if ref_var in sc_terms:
            print(f"There are sc_qss terms in {ref_var}...")
            print(sc_terms[ref_var])

            for scqss in sc_terms[ref_var]:
                print(scqss)
                term1 = sme.diff(sym_smp, scqss)
                scqssnum = self.syms_to_specnum(scqss)
                term2 = self.chain_diff(
                    species_info, self.sc_qss_smp[scqssnum], ref_var
                )

                dsym_dref += term1 * term2

                # print(species_info.dict_qssdepend_sc[species_info.qssa_species_list[scqssnum]])
                # print(species_info.dict_qssdepend_scqss[species_info.qssa_species_list[scqssnum]])
                # exit()
                # term2 = sme.diff(name, ref_var)
                # print(term2)
                # dsym_dref += sme.diff(sym_smp, name) * sme.diff(name, ref_var)

        return dsym_dref

    def compute_dscqss_dsc(self, scqss_idx, sc_idx, species_info):

        # Compute end of chain rule sc_qss derivatives
        self.compute_scqss_stopping(sc_idx, species_info)

        # print(self.dscqssdsc_stop)
        # print(len(self.dscqssdsc_stop))

        # exit()
        dscqss_dsc = 0
        return dscqss_dsc

    def compute_scqss_stopping(self, sc_idx, species_info):
        """Routine that computes the end terms in the sc_qss chain rule."""

        # Fill dscqssdsc_stop if not filled for sc yet
        if not self.dscqssdsc_stop["info"] == f"sc[{sc_idx}]":
            for stp in species_info.sc_qss_chain_stop:
                print(stp)
                print(species_info.dict_qss_species[stp])
                # if there is no dependence, then the derivative is zero
                if not species_info.dict_qssdepend_sc[stp]:
                    print(f"no dependence of {stp} on sc[{sc_idx}]")
                    tmp_diff = 0
                else:

                    times = time.time()
                    tmp_diff = sme.diff(
                        self.sc_qss_smp[species_info.dict_qss_species[stp]],
                        # smp.symbols(f"sc[{sc_idx}]"),
                        self.sc_smp[sc_idx],
                    )
                    print(
                        f"Time to do derivative for {stp} = {time.time()-times}"
                    )
                self.dscqssdsc_stop[
                    f"sc_qss[{species_info.dict_qss_species[stp]}]"
                ] = tmp_diff

            self.dscqssdsc_stop["info"] = f"sc[{sc_idx}]"
        else:
            # already filled for sc_idx...do nothing...
            print(f"dscqssdsc_stop already filled for {sc_idx}.")
