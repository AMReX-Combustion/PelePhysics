"""Generate C++ files for a mechanism."""
import pathlib
import shutil
import subprocess as spr

import numpy as np

import ceptr.ck as cck
import ceptr.formatter as cf
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
        jacobian=True,
        qss_format_input=None,
        qss_symbolic_jacobian=False,
    ):
        self.mechanism = mechanism

        self.jacobian = jacobian

        # Symbolic computations
        self.qss_symbolic_jacobian = qss_symbolic_jacobian

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
            qss_format_input,
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
                tempsp = csi.SpeciesDb(
                    id, sorted_idx, species.name, weight, species.charge
                )
                self.species_info.all_species.append(tempsp)
                self.species_info.nonqssa_species.append(tempsp)
                self.species_info.all_species_list.append(species.name)
                self.species_info.all_species_formatted_list.append(
                    cf.format_species(species.name)
                )
                self.species_info.nonqssa_species_list.append(species.name)
                self.species_info.nonqssa_species_formatted_list.append(
                    cf.format_species(species.name)
                )
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
                tempsp = csi.SpeciesDb(
                    id, sorted_idx, species.name, weight, species.charge
                )
                self.species_info.all_species.append(tempsp)
                self.species_info.qssa_species.append(tempsp)
                self.species_info.all_species_list.append(species.name)
                self.species_info.all_species_formatted_list.append(
                    cf.format_species(species.name)
                )
                self.species_info.qssa_species_list.append(species.name)
                self.species_info.qssa_species_formatted_list.append(
                    cf.format_species(species.name)
                )
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
        if self.species_info.n_qssa_species > 0:
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

            self.species_info.create_dicts()
            if self.species_info.n_qssa_species > 0:
                helper_names_to_print = []
                intermediate_names_to_print = []

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

                # self.species_info.create_dicts()
                self.species_info.identify_qss_dependencies(self.syms)
                self.species_info.identify_nonqss_dependencies(self.syms)
                self.species_info.make_scqss_dataframe()
                self.species_info.make_sc_dataframe()

                print(self.species_info.scqss_df)
                print(self.species_info.sc_df)

                # prod rate related
                cp.production_rate(
                    hdr,
                    self.mechanism,
                    self.species_info,
                    self.reaction_info,
                    self.syms,
                )
                cp.production_rate_light(
                    hdr,
                    self.mechanism,
                    self.species_info,
                    self.reaction_info,
                )

                # if self.species_info.n_qssa_species > 0:
                self.species_info.identify_wdot_dependencies(self.syms)
                self.species_info.make_wdot_dataframe()
                print(self.species_info.wdot_df)

                # Evaluate the dscqss_dscqss values for later
                self.syms.compute_dscqss_dscqss(species_info=self.species_info)

                # Evaluate the dscqss_dsc values for later
                self.syms.compute_dscqss_dsc(species_info=self.species_info)

                # # Evaluate the dwdot_dscqss values for later
                self.syms.compute_dwdot_dscqss(species_info=self.species_info)

                # # Evaluate the dwdot_dsc values for later
                self.syms.compute_dwdot_dsc(species_info=self.species_info)

                cck.ckwc(hdr, self.mechanism, self.species_info)
                cck.ckwyp(hdr, self.mechanism, self.species_info)
                cck.ckwxp(hdr, self.mechanism, self.species_info)
                cck.ckwyr(hdr, self.mechanism, self.species_info)
                cck.ckwxr(hdr, self.mechanism, self.species_info)
                cck.ckchrg(hdr, self)
                cck.ckchrgmass(hdr, self.species_info)
                cth.dthermodtemp(hdr, self.mechanism, self.species_info)

                # Approx analytical jacobian
                cj.ajac(
                    hdr,
                    self.mechanism,
                    self.species_info,
                    self.reaction_info,
                    jacobian=self.jacobian,
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
                # Analytical jacobian on GPU -- not used on CPU, define in mechanism.cpp
                if self.qss_symbolic_jacobian:
                    cj.ajac_symbolic(
                        hdr,
                        self.mechanism,
                        self.species_info,
                        self.reaction_info,
                        jacobian=self.jacobian,
                        syms=self.syms,
                    )
                else:
                    cj.ajac(
                        hdr,
                        self.mechanism,
                        self.species_info,
                        self.reaction_info,
                        jacobian=self.jacobian,
                    )
                cj.dproduction_rate(
                    hdr, self.mechanism, self.species_info, self.reaction_info
                )

            else:
                cp.production_rate(
                    hdr,
                    self.mechanism,
                    self.species_info,
                    self.reaction_info,
                    self.syms,
                )
                cck.ckwc(hdr, self.mechanism, self.species_info)
                cck.ckwyp(hdr, self.mechanism, self.species_info)
                cck.ckwxp(hdr, self.mechanism, self.species_info)
                cck.ckwyr(hdr, self.mechanism, self.species_info)
                cck.ckwxr(hdr, self.mechanism, self.species_info)
                cck.ckchrg(hdr, self)
                cck.ckchrgmass(hdr, self.species_info)
                cth.dthermodtemp(hdr, self.mechanism, self.species_info)
                # Approx analytical jacobian
                cj.ajac(
                    hdr,
                    self.mechanism,
                    self.species_info,
                    self.reaction_info,
                    jacobian=self.jacobian,
                    precond=True,
                )
                cj.dproduction_rate(
                    hdr,
                    self.mechanism,
                    self.species_info,
                    self.reaction_info,
                    precond=True,
                )
                # Analytical jacobian on GPU -- not used on CPU, define in mechanism.cpp
                cj.ajac(
                    hdr,
                    self.mechanism,
                    self.species_info,
                    self.reaction_info,
                    jacobian=self.jacobian,
                )
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
        if shutil.which(clexec) is not None:
            spr.run([clexec, "-i", self.hdrname])
            spr.run([clexec, "-i", self.cppname])
        else:
            print(
                "Clang-format not found. C++ files will be hard to parse by a"
                " human."
            )

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
                fstream, f"awt[{idx}] = {aw:f}; " + cw.comment(f"{elem}")
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
            text = f"imw_new[{i}] = {1.0 / species.weight:.16f};"
            cw.writer(fstream, text + cw.comment(f"{species.name}"))
        cw.writer(fstream, "}")
        cw.writer(fstream)

        cw.writer(fstream, cw.comment(" molecular weights "))
        cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
        cw.writer(fstream, "void get_mw(amrex::Real *mw_new){")
        for i in range(0, self.species_info.n_species):
            species = self.species_info.nonqssa_species[i]
            text = f"mw_new[{i}] = {species.weight:f};"
            cw.writer(fstream, text + cw.comment(f"{species.name}"))
        cw.writer(fstream, "}")

        cw.writer(fstream)
        cw.writer(fstream, cw.comment(" inverse molecular weights "))
        cw.writer(
            fstream,
            "AMREX_GPU_CONSTANT const amrex::Real "
            f"global_imw[{self.species_info.n_species}]={{",
        )
        for i in range(0, self.species_info.n_species):
            species = self.species_info.nonqssa_species[i]
            text = f"{1.0 / species.weight:.16f},"
            cw.writer(fstream, text + cw.comment(f"{species.name}"))
        cw.writer(fstream, "};")
        cw.writer(
            fstream,
            "const amrex::Real "
            f"h_global_imw[{self.species_info.n_species}]={{",
        )
        for i in range(0, self.species_info.n_species):
            species = self.species_info.nonqssa_species[i]
            text = f"{1.0 / species.weight:.16f},"
            cw.writer(fstream, text + cw.comment(f"{species.name}"))
        cw.writer(fstream, "};")
        cw.writer(fstream)

        cw.writer(fstream, cw.comment(" molecular weights "))
        cw.writer(
            fstream,
            "AMREX_GPU_CONSTANT const amrex::Real "
            f"global_mw[{self.species_info.n_species}]={{",
        )
        for i in range(0, self.species_info.n_species):
            species = self.species_info.nonqssa_species[i]
            text = f"{species.weight:f},"
            cw.writer(fstream, text + cw.comment(f"{species.name}"))
        cw.writer(fstream, "};")
        cw.writer(
            fstream,
            "const amrex::Real "
            f"h_global_mw[{self.species_info.n_species}]={{",
        )
        for i in range(0, self.species_info.n_species):
            species = self.species_info.nonqssa_species[i]
            text = f"{species.weight:f},"
            cw.writer(fstream, text + cw.comment(f"{species.name}"))
        cw.writer(fstream, "};")

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
            fstream, "void CKINU(const int i, int &nspec, int * ki, int * nu);"
        )
        cw.writer(
            fstream,
            "void CKKFKR(const amrex::Real P, const amrex::Real T,"
            + "const amrex::Real * x, amrex::Real *  q_f, amrex::Real *  q_r);",
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
            (
                "void SPARSITY_INFO_SYST(int * nJdata, const int * consP, int"
                " NCELLS);"
            ),
        )
        cw.writer(
            fstream,
            (
                "void SPARSITY_INFO_SYST_SIMPLIFIED(int * nJdata, const int *"
                " consP);"
            ),
        )
        cw.writer(
            fstream,
            (
                "void SPARSITY_PREPROC_CSC(int * rowVals, int * colPtrs, const"
                " int * consP, int NCELLS);"
            ),
        )
        cw.writer(
            fstream,
            (
                "void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const"
                " int * consP, int NCELLS, int base);"
            ),
        )
        cw.writer(
            fstream,
            (
                "void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtrs,"
                " const int * consP, int NCELLS, int base);"
            ),
        )
        cw.writer(
            fstream,
            (
                "void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int"
                " * colPtrs, int * indx, const int * consP);"
            ),
        )
        cw.writer(
            fstream,
            (
                "void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int"
                " * rowPtr, const int * consP, int base);"
            ),
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
            cw.writer(fstream, f"{self.mechanism.element_index(elem)}  {elem}")
            nb_elem += 1
        cw.writer(fstream, "*/")
        cw.writer(fstream)
        cw.writer(fstream, cw.comment("Species"))
        nb_ions = 0
        for species in self.species_info.nonqssa_species_list:
            s = cf.format_species(species)
            cw.writer(
                fstream,
                f"#define {s}_ID {self.species_info.ordered_idx_map[species]}",
            )
            if s[-1] == "n" or s[-1] == "p" or s == "E":
                nb_ions += 1
        cw.writer(fstream)
        cw.writer(fstream, f"#define NUM_ELEMENTS {nb_elem}")
        cw.writer(
            fstream, f"#define NUM_SPECIES {self.species_info.n_species}"
        )
        cw.writer(fstream, f"#define NUM_IONS {nb_ions}")
        cw.writer(
            fstream,
            f"#define NUM_REACTIONS {len(self.mechanism.reactions())}",
        )
        cw.writer(fstream)
        cw.writer(fstream, "#define NUM_FIT 4")
