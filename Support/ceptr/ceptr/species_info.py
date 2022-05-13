"""Species information."""
from collections import OrderedDict


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
        # non QSS
        # list of speciesDb for each non QSS spec
        self.nonqssa_species = []
        # list of non QSS species names
        self.nonqssa_species_list = []
        # number of non QSS species
        self.n_species = 0
        # QSS
        # list of speciesDb for each QSS spec
        self.qssa_species = []
        # list of QSS species names
        self.qssa_species_list = []
        # number of QSS species
        self.n_qssa_species = 0
        # all Species
        self.all_species = []
        self.all_species_list = []
        self.n_all_species = 0
        # Ordered dict for matching species to indices
        self.ordered_idx_map = OrderedDict()
        self.mech_idx_map = OrderedDict()

        self.low_temp = 0.0
        self.high_temp = 1000000.0

        # QSSA specific
        # sp-sp network
        self.qssa_ssnet = []
        # sp-reac network
        self.qssa_srnet = []
        # sp coupling network
        self.qssa_scnet = []
        # sp-sp network indices i of non zero elem
        self.qssa_ss_si = []
        # sp-sp network indices j of non zero elem
        self.qssa_ss_sj = []
        # sp-reac network indices i of non zero elem
        self.qssa_sr_si = []
        # sp-reac network indices j of non zero elem
        self.qssa_sr_sj = []
        # sp coupling network indices i of non zero elem
        self.qssa_sc_si = []
        # sp coupling network indices j of non zero elem
        self.qssa_sc_sj = []
