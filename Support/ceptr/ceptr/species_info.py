"""Species information."""
from collections import OrderedDict

import ceptr.qssa_info as cqi


class SpeciesDb:
    """Species database."""

    def __init__(self, idx, ordered_id, name, mwt):
        self.mech_idx = idx
        self.idx = ordered_id
        self.name = name
        self.weight = mwt

    def __str__(self):
        """Print members of SpeciesDb."""
        return (
            f"""name = {self.name}, """
            + f"""mech_idx = {self.mech_idx}, """
            + f"""idx = {self.idx}, """
            + f"""weight = {self.weight}"""
        )

    def __repr__(self):
        """Representation of SpeciesDb."""
        return f"""SpeciesDb({self.mech_idx}, {self.idx}, {self.name}, {self.weight})"""


class SpeciesInfo:
    """Information on species."""

    def __init__(self):
        # Species
        # non QSSA
        # list of speciesDb for each non QSSA spec
        self.nonqssa_species = []
        # list of non QSSA species names
        self.nonqssa_species_list = []
        # number of non QSSA species
        self.n_species = 0

        # all Species
        self.all_species = []
        self.all_species_list = []
        self.n_all_species = 0
        # Ordered dict for matching species to indices
        self.ordered_idx_map = OrderedDict()
        self.mech_idx_map = OrderedDict()

        self.low_temp = 0.0
        self.high_temp = 1000000.0

        # list of speciesDb for each QSSA spec
        self.qssa_species = []
        # list of QSSA species names
        self.qssa_species_list = []
        # number of QSSA species
        self.n_qssa_species = 0
        self.qssa_info = cqi.QSSAInfo()

    def create_dicts(self):
        """Create species dicts for ease of use."""
        self.dict_species = {v: i for i, v in enumerate(self.all_species_list)}
        self.dict_qss_species = {
            v: i for i, v in enumerate(self.qssa_species_list)
        }
        self.dict_nonqss_species = {
            v: i for i, v in enumerate(self.nonqssa_species_list)
        }

    def identify_qss_dependencies(self, syms):
        """Identify QSS species dependencies from syms object."""
        self.dict_qssdepend_scqss = {}
        self.dict_qssdepend_sc = {}
        self.dict_qssdepend_gRTqss = {}
        self.dict_qssdepend_gRT = {}
        self.dict_qssdepend_kf = {}
        self.dict_qssdepend_kr = {}
        for symbol in self.dict_qss_species:
            free_symb = syms.sc_qss_smp[
                self.dict_qss_species[symbol]
            ].free_symbols
            qss_symb = []
            sc_symb = []
            gRTqss_symb = []
            gRT_symb = []
            kf_symb = []
            kr_symb = []
            for ss in free_symb:
                if "sc_qss" in str(ss):
                    qss_symb.append(ss)
                elif "sc" in str(ss):
                    sc_symb.append(ss)
                elif "g_RT_qss" in str(ss):
                    gRTqss_symb.append(ss)
                elif "g_RT" in str(ss):
                    gRT_symb.append(ss)
                elif "kf" in str(ss):
                    kf_symb.append(ss)
                elif "kr" in str(ss):
                    kr_symb.append(ss)

            self.dict_qssdepend_scqss[symbol] = qss_symb
            self.dict_qssdepend_sc[symbol] = sc_symb
            self.dict_qssdepend_gRTqss[symbol] = gRTqss_symb
            self.dict_qssdepend_gRT[symbol] = gRT_symb
            self.dict_qssdepend_kf[symbol] = kf_symb
            self.dict_qssdepend_kr[symbol] = kr_symb
