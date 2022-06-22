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
            # Create dict to hold intermediate chain rule dscqssdsc terms
            self.dscqssdsc_interm = {}
            # Create dict to hold scqss_sc terms
            self.dscqssdsc = {}

    # @profile
    def convert_to_cpp(self, sym_smp):
        """Convert sympy object to C code compatible string."""
        # Convert to ccode (to fix pow) and then string
        cppcode = sme.ccode(sym_smp)
        cpp_str = str(cppcode)

        return cpp_str

    # @profile
    def write_array_to_cpp(
        self, list_smp, array_str, cw, fstream, indexList=None
    ):
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
            # Debugging prints in mechanism.H...
            # cw.writer(
            #     fstream,
            #     """std::cout << "%s = " << %s << std::endl;""" % (left_cse, left_cse),
            # )

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
            if indexList is None:
                cw.writer(
                    fstream,
                    "%s[%s] = %s;"
                    % (
                        array_str,
                        str(i),
                        cpp_str,
                    ),
                )
            else:
                cw.writer(
                    fstream,
                    "%s[%s] = %s;"
                    % (
                        array_str,
                        str(indexList[i]),
                        cpp_str,
                    ),
                )
            timee = time.time()
            print(
                "Printed expr for entry %d (time = %.3g s)"
                % (i, timee - times)
            )

    # @profile
    def syms_to_specnum(self, sym_smp):
        """Extracts number from syms string"""
        num = re.findall(r"\[(.*?)\]", str(sym_smp))
        return int(num[0])

    def compute_spec_dependencies(self, sym_smp):
        """Routine to compute the sc and sc_qss dependencies."""

        free_symb = sym_smp.free_symbols
        sc_depend = []
        scqss_depend = []
        for symbol in free_symb:
            if "sc_qss" in str(symbol):
                scqss_depend.append(symbol)
            elif "sc" in str(symbol):
                sc_depend.append(symbol)
            else:
                pass

        return sc_depend, scqss_depend

    # @profile
    def compute_dwdot_dsc(self, wdot_idx, sc_idx, species_info):
        """Routine to compute dwdot[x]/dsc[y]."""

        # First compute the dependencies of wdot
        print(f"Computing dependecies for wdot[{wdot_idx}]...")
        sc_depend, scqss_depend = self.compute_spec_dependencies(
            self.wdot_smp[wdot_idx]
        )

        # First compute dwdot/dsc[sc_idx]
        print(f"Computing dwdot[{wdot_idx}]/dsc[{sc_idx}]...")
        dwdotdsc = sme.diff(
            self.wdot_smp[wdot_idx], smp.symbols(f"sc[{sc_idx}]")
        )

        # Now compute dwdot/dscqss and the chain
        for scqss in scqss_depend:
            scqssnum = self.syms_to_specnum(scqss)
            print(f"Computing dwdot[{wdot_idx}]/dsc_qss[{scqssnum}]...")
            dwdotdscqss = sme.diff(
                self.wdot_smp[wdot_idx], smp.symbols(f"sc_qss[{scqssnum}]")
            )
            dscqss_dsc = self.compute_dscqss_dsc(
                scqssnum, sc_idx, species_info
            )
            dwdotdsc += dwdotdscqss * dscqss_dsc

        return dwdotdsc

    # @profile
    def compute_dscqss_dsc(self, scqss_idx, sc_idx, species_info):
        """Routine to compute dsc_qss[x]/dsc[y]."""

        if (scqss_idx, sc_idx) in self.dscqssdsc:
            # We have computed dscqssdsc before
            print(
                f"Returning dsc_qss[{scqss_idx}]/dsc[{sc_idx}] from memory..."
            )
            dscqss_dsc = self.dscqssdsc[(scqss_idx, sc_idx)]
        else:
            print(
                f"Compute stopping chain terms for scqss[{scqss_idx}]/sc[{sc_idx}]..."
            )
            self.compute_scqss_stopping(sc_idx, species_info)

            print(scqss_idx)
            debug_chain = f"dsc_qss[{scqss_idx}]/dsc[{sc_idx}] = "
            dscqss_dsc, debug_chain_out = self.chain_scqss_sc(
                scqss_idx, sc_idx, species_info
            )

            debug_chain += debug_chain_out

            # Store dscqss_dsc for later use...
            print(
                f"Storing dsc_qss[{scqss_idx}]/dsc[{sc_idx}] for later use..."
            )
            self.dscqssdsc[(scqss_idx, sc_idx)] = dscqss_dsc

            # print(debug_chain)

        return dscqss_dsc

    # @profile
    def chain_scqss_sc(self, scqss_idx, sc_idx, species_info):
        """Routine to compute chain rule scqss dependence recursively."""

        # Find the length of the sc_qss dependency list
        deplen = len(
            species_info.dict_qssdepend_scqss[
                species_info.qssa_species_list[scqss_idx]
            ]
        )
        print(f"deplen = {deplen}")
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

            if (scqss_idx, scqssnum) in self.dscqssdsc_interm:
                # We have computed that chain rule term before
                print(
                    f"Returning dsc_qss[{scqss_idx}]/dsc_qss[{scqssnum}] from memory..."
                )
                chain_vec[loop_idx] = self.dscqssdsc_interm[
                    (scqss_idx, scqssnum)
                ]
            else:
                print(f"Computing dsc_qss[{scqss_idx}]/dsc_qss[{scqssnum}]...")
                chain_vec[loop_idx] = sme.diff(
                    self.sc_qss_smp[scqss_idx],
                    smp.symbols(f"sc_qss[{scqssnum}]"),
                )
                self.dscqssdsc_interm[(scqss_idx, scqssnum)] = chain_vec[
                    loop_idx
                ]

            chain_vec_idx, chain_vec_debug_idx = self.chain_scqss_sc(
                scqssnum, sc_idx, species_info
            )
            # Multiply the result of the returned vectors by the current index
            chain_vec[loop_idx] *= chain_vec_idx
            chain_vec_debug[
                loop_idx
            ] = f" {chain_vec_debug[loop_idx]} * ({chain_vec_debug_idx}) "

            # self.dscqssdsc_interm[(scqss_idx, scqssnum)] = chain_vec[
            #     loop_idx
            # ]

        if deplen == 0:
            # If there are no dependencies, just return the end derivative
            if not self.dscqssdsc_stop["info"] == f"sc[{sc_idx}]":
                print(
                    f"dscqssdsc_stop should already be stored for {sc_idx}...!!!"
                )
                exit()
            else:
                print(
                    f"Returning pre-computed dscqssdsc_stop for sc_qss[{scqss_idx}]/sc[{sc_idx}]..."
                )
                chain_vec_out = self.dscqssdsc_stop[f"sc_qss[{scqss_idx}]"]
                chain_vec_debug_out = f"dsc_qss[{scqss_idx}]/dsc[{sc_idx}]"
        else:
            # sum up all terms for chain_vec
            chain_vec_out = 0
            for expr in chain_vec:
                chain_vec_out += expr
            #  Add in the symbolic derivative w.r.t. the initial sc term
            chain_vec_out += sme.diff(
                self.sc_qss_smp[scqss_idx], smp.symbols(f"sc[{sc_idx}]")
            )
            # sum up all terms for chain_vec_debug string
            chain_vec_debug_out = ""
            for item in chain_vec_debug:
                chain_vec_debug_out += f" + {item} "
            #  Add in the symbolic derivative w.r.t. the initial sc term
            chain_vec_debug_out += f" + dsc_qss[{scqss_idx}]/sc[{sc_idx}]"
        # Return both the computed sympy expressions and the debug string
        return chain_vec_out, chain_vec_debug_out

    # @profile
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
                    )
                    print(
                        f"Time to do derivative for scqss[{species_info.dict_qss_species[stp]}]/dsc[{sc_idx}] = {time.time()-times}"
                    )
                self.dscqssdsc_stop[
                    f"sc_qss[{species_info.dict_qss_species[stp]}]"
                ] = tmp_diff

            self.dscqssdsc_stop["info"] = f"sc[{sc_idx}]"
        else:
            # already filled for sc_idx...do nothing...
            print(f"dscqssdsc_stop already filled for {sc_idx}.")
            pass
