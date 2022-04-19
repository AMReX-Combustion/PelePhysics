from collections import OrderedDict


class SpeciesDb:
    """Species database."""

    def __init__(self, idx, ordered_id, name, mwt):
        self.mech_idx = idx
        self.idx = ordered_id
        self.name = name
        self.weight = mwt


class SpeciesInfo:
    """Information on species"""

    def __init__(self):
        # Species
        # non QSS
        # list of speciesDb for each non QSS spec
        self.nonqss_species = []
        # list of non QSS species names
        self.nonqss_species_list = []
        # number of non QSS species
        self.nSpecies = 0
        # QSS
        # list of speciesDb for each QSS spec
        self.qss_species = []
        # list of QSS species names
        self.qss_species_list = []
        # number of QSS species
        self.nQSSspecies = 0
        # all Species
        self.all_species = []
        self.all_species_list = []
        self.nAllspecies = 0
        # Ordered dict for matching species to indices
        self.ordered_idx_map = OrderedDict()
        self.mech_idx_map = OrderedDict()

        self.lowT = 0.0
        self.highT = 1000000.0
