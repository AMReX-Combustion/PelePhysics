"""Species information."""
from collections import OrderedDict

import ceptr.qssa_info as cqi


class SpeciesDb:
    """Species database."""

    def __init__(self, idx, ordered_id, name, mwt, chrg):
        self.mech_idx = idx
        self.idx = ordered_id
        self.name = name
        self.weight = mwt
        self.charge = chrg

    def __str__(self):
        """Print members of SpeciesDb."""
        return (
            f"""name = {self.name}, """
            + f"""mech_idx = {self.mech_idx}, """
            + f"""idx = {self.idx}, """
            + f"""weight = {self.weight}"""
            + f"""charge = {self.charge}"""
        )

    def __repr__(self):
        """Representation of SpeciesDb."""
        return f"""SpeciesDb({self.mech_idx}, {self.idx}, {self.name},
                             {self.weight}, {self.charge})"""


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
