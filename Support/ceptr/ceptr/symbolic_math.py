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
        cppcode = sme.ccode(sym_smp)
        cpp_str = str(cppcode)

        return cpp_str

    def write_array_to_cpp(self, list_smp, array_str, cw, fstream):
        """Convert sympy array to C code compatible string."""
        n = len(list_smp)

        # Write common expressions
        times = time.time()
        array_cse = sme.cse(list_smp)
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

    def compute_dscqss_dsc(self, scqss_idx, sc_idx, species_info):

        # Compute end of chain rule sc_qss derivatives
        self.compute_scqss_stopping(sc_idx, species_info)

        print(scqss_idx)
        debug_chain = f"dsc_qss[{scqss_idx}]/dsc[{sc_idx}] = "
        dscqss_dsc, debug_chain_out = self.chain_scqss(
            scqss_idx, sc_idx, species_info
        )

        debug_chain += debug_chain_out

        # print(dscqss_dsc)
        # print(debug_chain)
        # print(dscqss_dsc.free_symbols)

        return dscqss_dsc

    def chain_scqss(self, scqss_idx, sc_idx, species_info):

        # Find the length of the sc_qss dependency list
        deplen = len(
            species_info.dict_qssdepend_scqss[
                species_info.qssa_species_list[scqss_idx]
            ]
        )
        # Initialize vectors to store the recursive expressions
        chain_vec = [None] * deplen
        chain_vec_debug = [None] * deplen

        # Loop over the qss dependencies on scqss for the given scqss_idx
        for loop_idx, scqss_depend in enumerate(
            species_info.dict_qssdepend_scqss[
                species_info.qssa_species_list[scqss_idx]
            ]
        ):
            # Get the number of the scqss number from the syms obj
            scqssnum = self.syms_to_specnum(scqss_depend)
            # Make the debug string to ensure terms are properly computed
            chain_vec_debug[
                loop_idx
            ] = f"dsc_qss[{scqss_idx}]/dsc_qss[{scqssnum}]"
            # Compute the dsc_qss[scqss_idx]/dsc_qss[scqssnum] derivative
            print(f"Computing dsc_qss[{scqss_idx}]/dsc_qss[{scqssnum}]...")
            chain_vec[loop_idx] = sme.diff(
                self.sc_qss_smp[scqss_idx], smp.symbols(f"sc_qss[{scqssnum}]")
            )
            chain_vec_idx, chain_vec_debug_idx = self.chain_scqss(
                scqssnum, sc_idx, species_info
            )
            # Multiply the result of the returned vectors by the current index
            chain_vec[loop_idx] *= chain_vec_idx
            chain_vec_debug[
                loop_idx
            ] = f" {chain_vec_debug[loop_idx]} * ({chain_vec_debug_idx}) "

        if deplen == 0:
            # If there are no dependencies, just return the end derivative
            if not self.dscqssdsc_stop["info"] == f"sc[{sc_idx}]":
                print(
                    f"dscqssdsc_stop should already be stored for {sc_idx}...!!!"
                )
                exit()
            else:
                print(
                    f"Returning pre-computed dscqssdsc_stop for sc_qss[{scqss_idx}]"
                )
                chain_vec_out = self.dscqssdsc_stop[f"sc_qss[{scqss_idx}]"]
                chain_vec_debug_out = f"dsc_qss[{scqss_idx}]/dsc[{sc_idx}]"
        else:
            # sum up all terms for chain_vec
            chain_vec_out = 0
            for expr in chain_vec:
                chain_vec_out += expr
            # sum up all terms for chain_vec_debug string
            chain_vec_debug_out = ""
            for item in chain_vec_debug:
                chain_vec_debug_out += f" + {item} "
        # Return both the computed sympy expressions and the debug string
        return chain_vec_out, chain_vec_debug_out

    def compute_scqss_stopping(self, sc_idx, species_info):
        """Routine that computes the end terms in the sc_qss chain rule."""

        # Fill dscqssdsc_stop if not filled for sc yet
        if not self.dscqssdsc_stop["info"] == f"sc[{sc_idx}]":
            for stp in species_info.sc_qss_chain_stop:
                # if there is no dependence, then the derivative is zero
                if not species_info.dict_qssdepend_sc[stp]:
                    print(f"no dependence of {stp} on sc[{sc_idx}]")
                    tmp_diff = 0
                else:

                    times = time.time()
                    tmp_diff = sme.diff(
                        self.sc_qss_smp[species_info.dict_qss_species[stp]],
                        smp.symbols(f"sc[{sc_idx}]"),
                        # self.sc_smp[sc_idx],
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
            # print(f"dscqssdsc_stop already filled for {sc_idx}.")
            pass
