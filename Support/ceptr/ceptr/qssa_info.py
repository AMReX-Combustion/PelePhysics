"""QSSA information."""
from collections import OrderedDict


class QSSAInfo:
    """Information on QSSA."""

    def __init__(self):

        # sp-sp network
        self.ssnet = []
        # sp-reac network
        self.srnet = []
        # sp coupling network
        self.scnet = []
        # sp-sp network indices i of non zero elem
        self.ss_si = []
        # sp-sp network indices j of non zero elem
        self.ss_sj = []
        # sp-reac network indices i of non zero elem
        self.sr_si = []
        # sp-reac network indices j of non zero elem
        self.sr_sj = []
        # sp coupling network indices i of non zero elem
        self.sc_si = []
        # sp coupling network indices j of non zero elem
        self.sc_sj = []

        self.needs = OrderedDict()
        self.needs_count = OrderedDict()
        self.needs_running = OrderedDict()
        self.needs_count_running = OrderedDict()
        self.is_needed = OrderedDict()
        self.is_needed_count = OrderedDict()
        self.is_needed_running = OrderedDict()
        self.is_needed_count_running = OrderedDict()

        self.group = OrderedDict()
        self.discovery_order = 0
        self.potential_group = []
        self.lowest_link = OrderedDict()
        self.all_groups = OrderedDict()

        self.decouple_index = OrderedDict()
        self.decouple_count = 0

        self.rhs = OrderedDict()
        self.coeff = OrderedDict()
        self.qssa_coeff = OrderedDict()

        self.list_of_intermediate_helpers = []
