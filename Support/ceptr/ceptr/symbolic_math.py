"""Symbolic math for symbolic differentiation."""
import re
import time
from collections import OrderedDict

import numpy as np
import sympy as smp


class SymbolicMath:
    """Symbols to carry throughout operations."""

    def __init__(self, species_info, reaction_info):

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

    def syms_to_specnum(self, smp):
        """Extracts number from syms string"""
        num = re.findall(r"\[(.*?)\]", str(smp))
        return int(num[0])
