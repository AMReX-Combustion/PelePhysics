"""Generate C++ files for a mechanism."""
import pathlib
import shutil
import subprocess as spr
import sys

import ceptr.ck as cck
import ceptr.gjs as cgjs
import ceptr.jacobian as cj
import ceptr.production as cp
import ceptr.reaction_info as cri
import ceptr.sparsity as csp
import ceptr.species_info as csi
import ceptr.thermo as cth
import ceptr.transport as ctr
import ceptr.writer as cw


class Converter:
    """Convert Cantera mechanism to C++ files for Pele."""

    def __init__(self, mechanism):
        self.mechanism = mechanism
        self.mechpath = pathlib.Path(self.mechanism.source)
        self.hdrname = self.mechpath.parents[0] / "newmechanism.H"
        self.cppname = self.mechpath.parents[0] / "newmechanism.cpp"
        self.species_info = csi.SpeciesInfo()

        # Reactions
        # QSS
        self.qssReactions = []
        self.qfqr_co_idx_map = []
        self.nqssReactions = 0

        # QSS specific
        # sp-sp network
        self.QSS_SSnet = []
        # sp-reac network
        self.QSS_SRnet = []
        # sp coupling network
        self.QSS_SCnet = []
        # sp-sp network indices i of non zero elem
        self.QSS_SS_Si = []
        # sp-sp network indices j of non zero elem
        self.QSS_SS_Sj = []
        # sp-reac network indices i of non zero elem
        self.QSS_SR_Si = []
        # sp-reac network indices j of non zero elem
        self.QSS_SR_Rj = []
        # sp coupling network indices i of non zero elem
        self.QSS_SC_Si = []
        # sp coupling network indices j of non zero elem
        self.QSS_SC_Sj = []

        self.reacRemoveIDList = []
        try:
            f = open("reac_forward_to_remove", "r")
            lines = f.readlines()
            f.close()
            for line in lines:
                self.reacRemoveIDList.append(int(line))
        except FileNotFoundError:
            print("No forward reaction to remove")

        # List of intermediate helpers -- not optimal but can't be more clever rn
        self.list_of_intermediate_helpers = []

        self.setSpecies()
        # 0/ntroe/nsri/nlindem/nTB/nSimple/nWeird
        # 0/1    /2   /3      /4  /5      /6
        self.reaction_info = cri.sort_reactions(self.mechanism)

        # FIXME
        # # QSS  -- sort reactions/networks/check validity of QSSs
        # if self.nQSSspecies > 0:
        #     print("\n\n\n\n---------------------------------")
        #     print("+++++++++QSS INFO++++++++++++++++")
        #     print("---------------------------------")
        #     print("QSS species list =", self.qss_species_list)
        #     self._setQSSreactions(mechanism)
        #     self._getQSSnetworks(mechanism)  # sets up QSS subnetwork
        #     self._QSSvalidation(
        #         mechanism
        #     )  # Perform tests to ensure QSS species are good candidates
        #     self._QSSCoupling(
        #         mechanism
        #     )  # No quad coupling and fill SC network
        #     print("\n\n\n\n---------------------------------")
        #     print("+++++++++INIT NEEDS DICT+++++++++")
        #     print("---------------------------------")
        #     self._setQSSneeds(
        #         mechanism
        #     )  # Fill "need" dict (which species a species depends upon)
        #     self._setQSSisneeded(
        #         mechanism
        #     )  # Fill "is_needed" dict (which species needs that particular species)

    def setSpecies(self):
        """Set the species."""

        # Fill species counters
        self.species_info.nAllspecies = self.mechanism.n_species
        self.species_info.nQSSspecies = 0  # FIXME len(mechanism.qss_species())
        self.species_info.nSpecies = (
            self.species_info.nAllspecies - self.species_info.nQSSspecies
        )

        # FIXME QSS
        qss_list_tmp = []
        # # get the unsorted self.qss_species_list
        # for qss_sp in mechanism.qss_species():
        #     qss_list_tmp.append(qss_sp.symbol)

        # sort all species. First pass is for non QSS species
        # so we can put them at the beginning of the all species list
        sorted_idx = 0
        for id, species in enumerate(self.mechanism.species()):
            if species.name not in qss_list_tmp:
                weight = 0.0
                for elem, coef in species.composition.items():
                    aw = self.mechanism.atomic_weight(elem)
                    weight += coef * aw
                tempsp = csi.SpeciesDb(id, sorted_idx, species.name, weight)
                self.species_info.all_species.append(tempsp)
                self.species_info.nonqss_species.append(tempsp)
                self.species_info.all_species_list.append(species.name)
                self.species_info.nonqss_species_list.append(species.name)
                self.species_info.ordered_idx_map[species.name] = sorted_idx
                self.species_info.mech_idx_map[species.name] = id
                sorted_idx += 1

        # second pass through QSS species - put them at the end of the all spec list
        for id, species in enumerate(self.mechanism.species()):
            if species.name in qss_list_tmp:
                weight = 0.0
                for elem, coef in species.composition.items():
                    aw = self.mechanism.atomic_weight(elem)
                    weight += coef * aw
                tempsp = csi.SpeciesDb(id, sorted_idx, species.name, weight)
                self.species_info.all_species.append(tempsp)
                self.species_info.qss_species.append(tempsp)
                self.species_info.all_species_list.append(species.name)
                self.species_info.qss_species_list.append(species.name)
                self.species_info.ordered_idx_map[species.name] = sorted_idx
                self.species_info.mech_idx_map[species.name] = id
                sorted_idx += 1

        # FIXME
        # # Initialize QSS species-species, species-reaction, and species coupling networks
        # self.QSS_SSnet = np.zeros([self.nQSSspecies, self.nQSSspecies], "d")
        # self.QSS_SRnet = np.zeros(
        #     [self.nQSSspecies, self.mechanism.n_reactions)], "d"
        # )
        # self.QSS_SCnet = np.zeros([self.nQSSspecies, self.nQSSspecies], "d")

        print("FULL SPECIES LIST WITH TRANSPORTED FIRST AND QSS LAST: ")
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
            cw.writer(cpp, self.mechanism_includes())
            cri.rmap(cpp, self.mechanism, self.reaction_info)
            cri.get_rmap(cpp, self.mechanism)
            cck.ckinu(
                cpp, self.mechanism, self.species_info, self.reaction_info
            )
            self.atomicWeight(cpp)
            cck.ckawt(cpp, self.mechanism)
            cck.ckncf(cpp, self.mechanism, self.species_info)
            cck.cksyme_str(cpp, self.mechanism, self.species_info)
            cck.cksyms_str(cpp, self.mechanism, self.species_info)
            csp.sparsity(cpp, self.species_info)
            cw.writer(cpp, "#endif")

            # This is for the header file
            cw.writer(hdr, "#ifndef MECHANISM_H")
            cw.writer(hdr, "#define MECHANISM_H")
            self.print_mech_header(hdr)
            self.chem_file_decl(hdr)
            # Basic info
            cck.ckindx(hdr, self.mechanism, self.species_info)
            self.molecular_weights(hdr)
            cck.ckrp(hdr, self.mechanism, self.species_info)
            cth.thermo(hdr, self.mechanism, self.species_info)
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
            cck.T_given_ey(hdr)
            cck.T_given_hy(hdr)
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

            # QSS
            if self.species_info.nQSSspecies > 0:
                print("FIXME")
                sys.exit(1)
                # print("\n\n\n\n---------------------------------")
                # print("+++++++++GROUPS++++++++++++++++++")
                # print("---------------------------------")
                # self._getQSSgroups(mechanism)  # Figure out dependencies
                # print("\n\n\n\n---------------------------------")
                # print("+++++++++QSS SORTING+++++++++++++")
                # print("---------------------------------")
                # self._sortQSScomputation(
                #     mechanism
                # )  # Sort out order of group eval
                # print("\n\n\n\n---------------------------------")
                # print("+++++++++QSS EVAL++++++++++++++++")
                # print("---------------------------------")
                # self._sortQSSsolution_elements(
                #     mechanism
                # )  # Actually gauss-pivot the matrix to get algebraic expr
                # print("\n\n\n\n---------------------------------")
                # print("+++++++++QSS PRINTING++++++++++++")
                # print("---------------------------------")
                # self._QSScomponentFunctions(mechanism)

            # prod rate related
            cp.productionRate(
                hdr, self.mechanism, self.species_info, self.reaction_info
            )
            cck.ckwc(hdr, self.mechanism, self.species_info)
            cck.ckwyp(hdr, self.mechanism, self.species_info)
            cck.ckwxp(hdr, self.mechanism, self.species_info)
            cck.ckwyr(hdr, self.mechanism, self.species_info)
            cck.ckwxr(hdr, self.mechanism, self.species_info)
            cth.dthermodT(hdr, self.mechanism, self.species_info)
            # Approx analytical jacobian
            cj.ajacPrecond(
                hdr, self.mechanism, self.species_info, self.reaction_info
            )
            cj.DproductionRatePrecond(
                hdr, self.mechanism, self.species_info, self.reaction_info
            )
            # # Analytical jacobian on GPU -- not used on CPU, define in mechanism.cpp
            cj.ajac(hdr, self.mechanism, self.species_info, self.reaction_info)
            cj.DproductionRate(
                hdr, self.mechanism, self.species_info, self.reaction_info
            )
            # Transport
            cw.writer(hdr)
            ctr.transport(hdr, self.mechanism, self.species_info)
            ctr.getCriticalParameters(hdr, self.mechanism, self.species_info)
            # GS routines
            cgjs.emptygjs(hdr)
            cw.writer(hdr)
            cw.writer(hdr, "#endif")

    def mechanism_includes(self):
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

    def atomicWeight(self, fstream):
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
        cw.writer(fstream)
        cw.writer(fstream, cw.comment(" inverse molecular weights "))
        cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
        cw.writer(fstream, "void get_imw(amrex::Real *imw_new){")
        for i in range(0, self.species_info.nSpecies):
            species = self.species_info.nonqss_species[i]
            text = "imw_new[%d] = 1.0/%f;" % (i, species.weight)
            cw.writer(fstream, text + cw.comment("%s" % species.name))
        cw.writer(fstream, "}")
        cw.writer(fstream)

        cw.writer(fstream, cw.comment(" molecular weights "))
        cw.writer(fstream, "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE")
        cw.writer(fstream, "void get_mw(amrex::Real *mw_new){")
        for i in range(0, self.species_info.nSpecies):
            species = self.species_info.nonqss_species[i]
            text = "mw_new[%d] = %f;" % (i, species.weight)
            cw.writer(fstream, text + cw.comment("%s" % species.name))
        cw.writer(fstream, "}")

    def chem_file_decl(self, fstream):
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

    def print_mech_header(self, fstream):
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
        for species in self.species_info.nonqss_species_list:
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
            fstream, "#define NUM_SPECIES %d" % (self.species_info.nSpecies)
        )
        cw.writer(
            fstream,
            "#define NUM_REACTIONS %d" % (len(self.mechanism.reactions())),
        )
        cw.writer(fstream)
        cw.writer(fstream, "#define NUM_FIT 4")
