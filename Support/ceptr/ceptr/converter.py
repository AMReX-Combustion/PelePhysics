"""Generate C++ files for a mechanism."""
import pathlib
import re
import shutil
import subprocess as spr
import time

import numpy as np
import symengine as sme
import sympy as smp

import ceptr.ck as cck
import ceptr.gjs as cgjs
import ceptr.jacobian as cj
import ceptr.production as cp
import ceptr.qssa_converter as cqc
import ceptr.reaction_info as cri
import ceptr.sparsity as csp
import ceptr.species_info as csi
import ceptr.symbolic_math as csm
import ceptr.thermo as cth
import ceptr.transport as ctr
import ceptr.writer as cw


class Converter:
    """Convert Cantera mechanism to C++ files for Pele."""

    def __init__(
        self,
        mechanism,
        hformat,
        remove_1,
        remove_pow2,
        min_op_count,
        recursive_op_count,
        store_in_jacobian,
        round_decimals,
    ):
        self.mechanism = mechanism
        self.hformat = hformat
        self.remove_1 = remove_1
        self.remove_pow2 = remove_pow2
        self.min_op_count = min_op_count
        self.recursive_op_count = recursive_op_count
        self.store_in_jacobian = store_in_jacobian
        self.round_decimals = round_decimals
        self.mechpath = pathlib.Path(self.mechanism.source)
        self.rootname = "mechanism"
        self.hdrname = self.mechpath.parents[0] / f"{self.rootname}.H"
        self.cppname = self.mechpath.parents[0] / f"{self.rootname}.cpp"
        self.species_info = csi.SpeciesInfo()

        self.set_species()
        # 0/ntroe/nsri/nlindem/nTB/nSimple/nWeird
        # 0/1    /2   /3      /4  /5      /6
        self.reaction_info = cri.sort_reactions(self.mechanism)
        # QSS  -- sort reactions/networks/check validity of QSSs
        if self.species_info.n_qssa_species > 0:
            print("QSSA information")
            print("QSS species list =", self.species_info.qssa_species_list)
            cqc.set_qssa_reactions(
                self.mechanism, self.species_info, self.reaction_info
            )
            cqc.get_qssa_networks(
                self.mechanism, self.species_info, self.reaction_info
            )
            # sets up QSS subnetwork
            cqc.get_qssa_networks(
                self.mechanism, self.species_info, self.reaction_info
            )
            # Perform tests to ensure QSSA species are good candidates
            cqc.qssa_validation(
                self.mechanism, self.species_info, self.reaction_info
            )
            # No quad coupling and fill SC network
            cqc.qssa_coupling(
                self.mechanism, self.species_info, self.reaction_info
            )
            # Fill "need" dict (which species a species depends upon)
            print("QSSA initialization needs dictionary")
            cqc.set_qssa_needs(
                self.mechanism, self.species_info, self.reaction_info
            )
            # Fill "is_needed" dict (which species needs that particular species)
            cqc.set_qssa_isneeded(
                self.mechanism, self.species_info, self.reaction_info
            )
        # Initialize symbolic variables
        self.syms = csm.SymbolicMath(
            self.species_info,
            self.reaction_info,
            self.mechanism,
            self.hformat,
            self.remove_1,
            self.remove_pow2,
            self.min_op_count,
            self.recursive_op_count,
            self.store_in_jacobian,
            self.round_decimals,
        )

    def set_species(self):
        """Set the species."""
        # Fill species counters
        self.species_info.n_all_species = self.mechanism.n_species
        try:
            self.species_info.n_qssa_species = self.mechanism.input_data[
                "n_qssa_species"
            ]
        except KeyError:
            self.species_info.n_qssa_species = 0

        self.species_info.n_species = (
            self.species_info.n_all_species - self.species_info.n_qssa_species
        )

        # get the unsorted self.qssa_species_list
        qssa_list_tmp = []
        try:
            for qssa_sp in self.mechanism.input_data["qssa_species"]:
                qssa_list_tmp.append(qssa_sp)
        except KeyError:
            pass

        # sort all species. First pass is for non QSS species
        # so we can put them at the beginning of the all species list
        sorted_idx = 0
        for id, species in enumerate(self.mechanism.species()):
            if species.name not in qssa_list_tmp:
                weight = 0.0
                for elem, coef in species.composition.items():
                    aw = self.mechanism.atomic_weight(elem)
                    weight += coef * aw
                tempsp = csi.SpeciesDb(id, sorted_idx, species.name, weight)
                self.species_info.all_species.append(tempsp)
                self.species_info.nonqssa_species.append(tempsp)
                self.species_info.all_species_list.append(species.name)
                self.species_info.nonqssa_species_list.append(species.name)
                self.species_info.ordered_idx_map[species.name] = sorted_idx
                self.species_info.mech_idx_map[species.name] = id
                sorted_idx += 1

        # second pass through QSS species - put them at the end of the all spec list
        for id, species in enumerate(self.mechanism.species()):
            if species.name in qssa_list_tmp:
                weight = 0.0
                for elem, coef in species.composition.items():
                    aw = self.mechanism.atomic_weight(elem)
                    weight += coef * aw
                tempsp = csi.SpeciesDb(id, sorted_idx, species.name, weight)
                self.species_info.all_species.append(tempsp)
                self.species_info.qssa_species.append(tempsp)
                self.species_info.all_species_list.append(species.name)
                self.species_info.qssa_species_list.append(species.name)
                self.species_info.ordered_idx_map[species.name] = sorted_idx
                self.species_info.mech_idx_map[species.name] = id
                sorted_idx += 1

        # Initialize QSS species-species, species-reaction, and species coupling networks
        self.species_info.qssa_info.ssnet = np.zeros(
            [
                self.species_info.n_qssa_species,
                self.species_info.n_qssa_species,
            ],
            "d",
        )
        self.species_info.qssa_info.srnet = np.zeros(
            [self.species_info.n_qssa_species, self.mechanism.n_reactions], "d"
        )
        self.species_info.qssa_info.scnet = np.zeros(
            [
                self.species_info.n_qssa_species,
                self.species_info.n_qssa_species,
            ],
            "d",
        )

        print("Full species list with transported first and QSSA last:")
        for all_species in self.species_info.all_species:
            print(
                all_species.name,
                " ",
                all_species.idx,
                " ",
                all_species.mech_idx,
                " ",
                all_species.weight,
            )

    def writer(self):
        """Write out the C++ files."""
        with open(self.hdrname, "w") as hdr, open(self.cppname, "w") as cpp:
            # This is for the cpp file
            cw.writer(cpp, self.mechanism_cpp_includes())
            cri.rmap(cpp, self.mechanism, self.reaction_info)
            cri.get_rmap(cpp, self.mechanism)
            cck.ckinu(
                cpp, self.mechanism, self.species_info, self.reaction_info
            )
            cck.ckkfkr(cpp, self.mechanism, self.species_info)
            cp.progress_rate_fr(
                cpp, self.mechanism, self.species_info, self.reaction_info
            )
            self.atomic_weight(cpp)
            cck.ckawt(cpp, self.mechanism)
            cck.ckncf(cpp, self.mechanism, self.species_info)
            cck.cksyme_str(cpp, self.mechanism, self.species_info)
            cck.cksyms_str(cpp, self.mechanism, self.species_info)
            csp.sparsity(cpp, self.species_info)

            # This is for the header file
            cw.writer(hdr, "#ifndef MECHANISM_H")
            cw.writer(hdr, "#define MECHANISM_H")
            self.mechanism_header_includes(hdr)
            self.mechanism_cpp_declarations(hdr)
            # Basic info
            cck.ckindx(hdr, self.mechanism, self.species_info)
            self.molecular_weights(hdr)
            cck.ckrp(hdr, self.mechanism, self.species_info)
            cth.thermo(hdr, self.mechanism, self.species_info, self.syms)
            # mean quantities -- do not take QSS into account, sumX and Y = 1 without them
            cck.ckcpbl(hdr, self.mechanism, self.species_info)
            cck.ckcpbs(hdr, self.mechanism, self.species_info)
            cck.ckcvbl(hdr, self.mechanism, self.species_info)
            cck.ckcvbs(hdr, self.mechanism, self.species_info)
            cck.ckhbml(hdr, self.mechanism, self.species_info)
            cck.ckhbms(hdr, self.mechanism, self.species_info)
            cck.ckubml(hdr, self.mechanism, self.species_info)
            cck.ckubms(hdr, self.mechanism, self.species_info)
            cck.cksbml(hdr, self.mechanism, self.species_info)
            cck.cksbms(hdr, self.mechanism, self.species_info)
            cck.temp_given_ey(hdr)
            cck.temp_given_hy(hdr)
            cck.ckpx(hdr, self.mechanism, self.species_info)
            cck.ckpy(hdr, self.mechanism, self.species_info)
            cck.ckpc(hdr, self.mechanism, self.species_info)
            cck.ckrhox(hdr, self.mechanism, self.species_info)
            cck.ckrhoy(hdr, self.mechanism, self.species_info)
            cck.ckrhoc(hdr, self.mechanism, self.species_info)
            cck.ckwt(hdr, self.mechanism, self.species_info)
            cck.ckmmwy(hdr, self.mechanism, self.species_info)
            cck.ckmmwx(hdr, self.mechanism, self.species_info)
            cck.ckmmwc(hdr, self.mechanism, self.species_info)
            cck.ckcpor(hdr, self.mechanism, self.species_info)
            cck.ckhort(hdr, self.mechanism, self.species_info)
            cck.cksor(hdr, self.mechanism, self.species_info)
            # conversions
            cck.ckytx(hdr, self.mechanism, self.species_info)
            cck.ckytcp(hdr, self.mechanism, self.species_info)
            cck.ckytcr(hdr, self.mechanism, self.species_info)
            cck.ckxty(hdr, self.mechanism, self.species_info)
            cck.ckxtcp(hdr, self.mechanism, self.species_info)
            cck.ckxtcr(hdr, self.mechanism, self.species_info)
            cck.ckctx(hdr, self.mechanism, self.species_info)
            cck.ckcty(hdr, self.mechanism, self.species_info)
            # species quantities
            # MOL
            cck.ckcvml(hdr, self.mechanism, self.species_info)
            cck.ckcpml(hdr, self.mechanism, self.species_info)
            cck.ckuml(hdr, self.mechanism, self.species_info)
            cck.ckhml(hdr, self.mechanism, self.species_info)
            # cck.ckgml(hdr, self.mechanism, self.species_info)
            # cck.ckaml(hdr, self.mechanism, self.species_info)
            cck.cksml(hdr, self.mechanism, self.species_info)
            # MASS
            cck.ckcvms(hdr, self.mechanism, self.species_info)
            cck.ckcpms(hdr, self.mechanism, self.species_info)
            cck.ckums(hdr, self.mechanism, self.species_info)
            cck.ckhms(hdr, self.mechanism, self.species_info)
            # cck.ckgms(hdr, self.mechanism, self.species_info)
            # cck.ckams(hdr, self.mechanism, self.species_info)
            cck.cksms(hdr, self.mechanism, self.species_info)

            if self.species_info.n_qssa_species > 0:

                helper_names_to_print = ["H_2"]
                intermediate_names_to_print = ["PXC5H11_rhs", "PXC7H15_rhs"]

                print("QSSA groups")
                # Figure out dependencies
                cqc.get_qssa_groups(
                    self.mechanism, self.species_info, self.reaction_info
                )
                print("QSSA sorting")
                # Sort out order of group evaluation
                cqc.sort_qssa_computation(
                    self.mechanism, self.species_info, self.reaction_info
                )
                # Invert QSSA print coeff and QSSA evaluation to see expressions
                #  in terms of qr and qf
                print("QSSA print coeff")
                cqc.qssa_coeff_functions(
                    hdr,
                    self.mechanism,
                    self.species_info,
                    self.reaction_info,
                    self.syms,
                )
                print("QSSA evaluation")
                # Actually gauss-pivot the matrix to get algebraic expr
                cqc.sort_qssa_solution_elements(
                    self.mechanism,
                    self.species_info,
                    self.reaction_info,
                    self.syms,
                )
                print("QSSA printing")
                cqc.qssa_component_functions(
                    hdr,
                    self.mechanism,
                    self.species_info,
                    self.reaction_info,
                    self.syms,
                    helper_names_to_print,
                    intermediate_names_to_print,
                )
                # print("Symbolic kf QSS print for debug")
                # cqc.qssa_kf_debug(
                #    hdr,
                #    self.mechanism,
                #    self.species_info,
                #    self.reaction_info,
                #    self.syms,
                # )
                # print("Symbolic thermo QSS print for debug")
                # cth.gibbsQSS_debug(
                #    hdr,
                #    self.mechanism,
                #    self.species_info,
                #    self.reaction_info,
                #    self.syms,
                # )
                # cth.speciesEnthalpyQSS_debug(
                #    hdr,
                #    self.mechanism,
                #    self.species_info,
                #    self.reaction_info,
                #    self.syms,
                # )
            #    print("Symbolic Sc qss print for debug")
            #    cqc.qssa_sc_qss_debug(
            #        hdr,
            #        self.mechanism,
            #        self.species_info,
            #        self.reaction_info,
            #        self.syms,
            #    )
            #    print("Symbolic qf qss print for debug")
            #    cqc.qssa_coeff_debug(
            #        hdr,
            #        self.mechanism,
            #        self.species_info,
            #        self.reaction_info,
            #        self.syms,
            #    )
            #    print("Symbolic qss terms print for debug")
            #    cqc.qssa_terms_debug(
            #        hdr,
            #        self.mechanism,
            #        self.species_info,
            #        self.reaction_info,
            #        self.syms,
            #        helper_names_to_print,
            #        intermediate_names_to_print,
            #    )

            # print("Symbolic thermo print for debug")
            # cth.gibbs_debug(
            #    hdr,
            #    self.mechanism,
            #    self.species_info,
            #    self.reaction_info,
            #    self.syms,
            # )
            # cth.speciesEnthalpy_debug(
            #    hdr,
            #    self.mechanism,
            #    self.species_info,
            #    self.reaction_info,
            #    self.syms,
            # )

            self.species_info.create_dicts()
            self.species_info.identify_qss_dependencies(self.syms)
            self.species_info.identify_nonqss_dependencies(self.syms)
            self.species_info.make_scqss_dataframe()
            self.species_info.make_sc_dataframe()

            print(self.species_info.scqss_df)
            print(self.species_info.sc_df)

            # prod rate related
            times = time.time()
            cp.production_rate(
                hdr,
                self.mechanism,
                self.species_info,
                self.reaction_info,
                self.syms,
            )
            print(f"Time to do production_rate = {time.time()-times}")
            times = time.time()
            cp.production_rate_light(
                hdr,
                self.mechanism,
                self.species_info,
                self.reaction_info,
            )
            print(f"Time to do production_rate light = {time.time()-times}")

            # print("Symbolic wdot print for debug")
            # cp.production_rate_debug(
            #    hdr,
            #    self.mechanism,
            #    self.species_info,
            #    self.reaction_info,
            #    self.syms,
            # )

            times = time.time()
            self.species_info.identify_wdot_dependencies(self.syms)
            self.species_info.make_wdot_dataframe()
            print(
                f"Time to identify wdot dependencies and make dataframe = {time.time()-times}"
            )
            print(self.species_info.wdot_df)

            # Evaluate the dscqss_dscqss values for later
            times = time.time()
            self.syms.compute_dscqss_dscqss(species_info=self.species_info)
            print(f"Time to do all dscqss_dscqss = {time.time()-times}")

            # Evaluate the dscqss_dsc values for later
            times = time.time()
            self.syms.compute_dscqss_dsc_fast(species_info=self.species_info)
            print(f"Time to do all the dscqss_dsc = {time.time()-times}")

            # # Evaluate the dwdot_dscqss values for later
            times = time.time()
            self.syms.compute_dwdot_dscqss_fast(species_info=self.species_info)
            print(f"Time to do all the dwdot_dscqss = {time.time()-times}")

            # # Evaluate the dwdot_dsc values for later
            times = time.time()
            self.syms.compute_dwdot_dsc_fast(species_info=self.species_info)
            print(f"Time to do all the dwdot_dsc = {time.time()-times}")

            cck.ckwc(hdr, self.mechanism, self.species_info)
            cck.ckwyp(hdr, self.mechanism, self.species_info)
            cck.ckwxp(hdr, self.mechanism, self.species_info)
            cck.ckwyr(hdr, self.mechanism, self.species_info)
            cck.ckwxr(hdr, self.mechanism, self.species_info)
            cth.dthermodtemp(hdr, self.mechanism, self.species_info)

            # print("Symbolic dscqss_dsc term print for debug")
            # cj.dscqss_dsc_debug(
            #     hdr,
            #     self.mechanism,
            #     self.species_info,
            #     self.reaction_info,
            #     self.syms,
            #     [
            #         dscqss0dsc0,
            #         dscqss1dsc0,
            #         dscqss2dsc0,
            #     ],
            #     [
            #         (self.species_info.n_species) * 0 + 0,
            #         (self.species_info.n_species) * 0 + 1,
            #         (self.species_info.n_species) * 0 + 2,
            #     ],
            # )

            # print("Symbolic dscqss_dsc term print for debug")
            # cj.dscqss_dsc_fast_debug(
            #     hdr,
            #     self.mechanism,
            #     self.species_info,
            #     self.reaction_info,
            #     self.syms,
            # )

            cj.ajac_term_fast_debug(
                hdr,
                self.mechanism,
                self.species_info,
                self.reaction_info,
                self.syms,
            )

            # Approx analytical jacobian
            cj.ajac(
                hdr,
                self.mechanism,
                self.species_info,
                self.reaction_info,
                precond=True,
                syms=self.syms,
            )
            cj.dproduction_rate(
                hdr,
                self.mechanism,
                self.species_info,
                self.reaction_info,
                precond=True,
            )
            # # Analytical jacobian on GPU -- not used on CPU, define in mechanism.cpp
            # cj.ajac(
            #     hdr,
            #     self.mechanism,
            #     self.species_info,
            #     self.reaction_info,
            #     syms=self.syms,
            # )
            cj.dproduction_rate(
                hdr, self.mechanism, self.species_info, self.reaction_info
            )
            # Transport
            cw.writer(hdr)
            ctr.transport(hdr, self.mechanism, self.species_info)
            ctr.critical_parameters(hdr, self.mechanism, self.species_info)
            # GS routines
            cgjs.emptygjs(hdr)
            cw.writer(hdr)
            cw.writer(hdr, "#endif")

    def mechanism_cpp_includes(self):
        """Write the mechanism cpp includes."""
        return '#include "mechanism.H"'

    def formatter(self):
        """Format with clang-format."""
        clexec = "clang-format"
        try:
            shutil.which(clexec)
        except shutil.Error:
            print("Clang-format not found")

        spr.run([clexec, "-i", self.hdrname])
        spr.run([clexec, "-i", self.cppname])

    def atomic_weight(self, fstream):
        """Write the atomic weight."""
        cw.writer(fstream)
        cw.writer(fstream, cw.comment("save atomic weights into array"))
        cw.writer(fstream, "void atomicWeight(amrex::Real *  awt)")
        cw.writer(fstream, "{")
        for elem in self.mechanism.element_names:
            idx = self.mechanism.element_index(elem)
            aw = self.mechanism.atomic_weight(elem)
            cw.writer(
                fstream, "awt[%d] = %f; " % (idx, aw) + cw.comment("%s" % elem)
            )
        cw.writer(fstream, "}")

    def molecular_weights(self, fstream):
        """Write the molecular weights."""
        cw.writer(fstream)
        cw.writer(fstream, cw.comment(" inverse molecular weights "))
        cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
        cw.writer(fstream, "void get_imw(amrex::Real *imw_new){")
        for i in range(0, self.species_info.n_species):
            species = self.species_info.nonqssa_species[i]
            text = "imw_new[%d] = 1.0/%f;" % (i, species.weight)
            cw.writer(fstream, text + cw.comment("%s" % species.name))
        cw.writer(fstream, "}")
        cw.writer(fstream)

        cw.writer(fstream, cw.comment(" molecular weights "))
        cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
        cw.writer(fstream, "void get_mw(amrex::Real *mw_new){")
        for i in range(0, self.species_info.n_species):
            species = self.species_info.nonqssa_species[i]
            text = "mw_new[%d] = %f;" % (i, species.weight)
            cw.writer(fstream, text + cw.comment("%s" % species.name))
        cw.writer(fstream, "}")

    def mechanism_cpp_declarations(self, fstream):
        """Write the chemistry function declarations."""
        cw.writer(fstream)
        cw.writer(
            fstream,
            cw.comment(
                " ALWAYS on CPU stuff -- can have different def depending on"
                " if we are CPU or GPU based. Defined in mechanism.cpp "
            ),
        )
        cw.writer(fstream, "void atomicWeight(amrex::Real *  awt);")
        cw.writer(fstream, cw.comment(" MISC "))
        cw.writer(fstream, "void CKAWT(amrex::Real *  awt);")
        cw.writer(fstream, "void CKNCF(int * ncf);")
        cw.writer(
            fstream, "void CKSYME_STR(amrex::Vector<std::string>& ename);"
        )
        cw.writer(
            fstream, "void CKSYMS_STR(amrex::Vector<std::string>& kname);"
        )
        cw.writer(fstream, "void GET_RMAP(int * _rmap);")
        cw.writer(
            fstream, "void CKINU(int * i, int * nspec, int * ki, int * nu);"
        )
        cw.writer(
            fstream,
            "void CKKFKR(amrex::Real *  P, amrex::Real *  T,"
            + "amrex::Real *  x, amrex::Real *  q_f, amrex::Real *  q_r);",
        )
        cw.writer(
            fstream,
            "void progressRateFR(amrex::Real *  q_f, amrex::Real *  q_r,"
            + "amrex::Real *  sc, amrex::Real T);",
        )
        cw.writer(fstream, cw.comment(" SPARSE INFORMATION "))
        cw.writer(
            fstream,
            "void SPARSITY_INFO(int * nJdata, const int * consP, int NCELLS);",
        )
        cw.writer(
            fstream,
            "void SPARSITY_INFO_SYST(int * nJdata, const int * consP, int"
            " NCELLS);",
        )
        cw.writer(
            fstream,
            "void SPARSITY_INFO_SYST_SIMPLIFIED(int * nJdata, const int *"
            " consP);",
        )
        cw.writer(
            fstream,
            "void SPARSITY_PREPROC_CSC(int * rowVals, int * colPtrs, const int"
            " * consP, int NCELLS);",
        )
        cw.writer(
            fstream,
            "void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int"
            " * consP, int NCELLS, int base);",
        )
        cw.writer(
            fstream,
            "void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtrs,"
            " const int * consP, int NCELLS, int base);",
        )
        cw.writer(
            fstream,
            "void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int *"
            " colPtrs, int * indx, const int * consP);",
        )
        cw.writer(
            fstream,
            "void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int *"
            " rowPtr, const int * consP, int base);",
        )

    def mechanism_header_includes(self, fstream):
        """Write the mechanism header includes."""
        cw.writer(fstream)
        cw.writer(fstream, "#include <AMReX_Gpu.H>")
        cw.writer(fstream, "#include <AMReX_REAL.H>")
        cw.writer(fstream)
        cw.writer(fstream, "/* Elements")
        nb_elem = 0
        for elem in self.mechanism.element_names:
            cw.writer(
                fstream, "%d  %s" % (self.mechanism.element_index(elem), elem)
            )
            nb_elem += 1
        cw.writer(fstream, "*/")
        cw.writer(fstream)
        cw.writer(fstream, cw.comment("Species"))
        for species in self.species_info.nonqssa_species_list:
            s = species.strip()
            # Ionic species
            if s[-1] == "-":
                s = s[:-1] + "n"
            if s[-1] == "+":
                s = s[:-1] + "p"
            # Excited species
            s = s.replace("*", "D")
            # Remove other characters not allowed in preprocessor defines
            s = s.replace("-", "").replace("(", "").replace(")", "")
            cw.writer(
                fstream,
                "#define %s_ID %d"
                % (s, self.species_info.ordered_idx_map[species]),
            )
        cw.writer(fstream)
        cw.writer(fstream, "#define NUM_ELEMENTS %d" % (nb_elem))
        cw.writer(
            fstream, "#define NUM_SPECIES %d" % (self.species_info.n_species)
        )
        cw.writer(
            fstream,
            "#define NUM_REACTIONS %d" % (len(self.mechanism.reactions())),
        )
        cw.writer(fstream)
        cw.writer(fstream, "#define NUM_FIT 4")
