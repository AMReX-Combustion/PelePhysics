"""Reaction information."""

from collections import OrderedDict

import ceptr.constants as cc
import ceptr.writer as cw


class ReactionInfo:
    """Information on reactions."""

    def __init__(self, mechanism, interface):
        """Initialize the reaction information."""
        self.n_reactions = mechanism.n_reactions
        self.rs = []
        self.rs_unsorted = mechanism.reactions()
        self.is_sorted = False
        self.index = [0]
        self.idxmap = OrderedDict()
        self.remove_id_list = []
        try:
            self.remove_id_list = mechanism.input_data["forward_to_remove_idx"]
        except KeyError:
            pass

        self.qssa_reactions = []
        self.n_qssa_reactions = 0
        self.qfqr_co_idx_map = []

        if not isinstance(interface, type(None)):
            self.rs_unsorted += interface.reactions()

def sort_reactions(mechanism, interface):
    """Sort reactions."""
    reaction_info = ReactionInfo(mechanism, interface)
    i = 0

    # troe
    for k, r in enumerate(reaction_info.rs_unsorted):
        if r not in reaction_info.rs:
            if r.reaction_type == "falloff-Troe":
                reaction_info.idxmap[k] = i
                reaction_info.rs.append(r)
                i += 1
    reaction_info.index.append(i)

    # sri
    for k, r in enumerate(reaction_info.rs_unsorted):
        if r not in reaction_info.rs:
            if r.reaction_type == "falloff-Sri":
                reaction_info.idxmap[k] = i
                reaction_info.rs.append(r)
                i += 1
    reaction_info.index.append(i)

    # lindemann
    for k, r in enumerate(reaction_info.rs_unsorted):
        if r not in reaction_info.rs:
            if r.reaction_type == "falloff-Lindemann":
                reaction_info.idxmap[k] = i
                reaction_info.rs.append(r)
                i += 1
    reaction_info.index.append(i)

    # three-body
    for k, r in enumerate(reaction_info.rs_unsorted):
        if r not in reaction_info.rs:
            if r.third_body is not None:
                reaction_info.idxmap[k] = i
                reaction_info.rs.append(r)
                i += 1
    reaction_info.index.append(i)

    # simplest case (Arrhenius reactions)
    for k, r in enumerate(reaction_info.rs_unsorted):
        if r not in reaction_info.rs:
            # skip all interface-Arrhenius and sticking-Arrhenius heterogeneous reactions
            # heterogeneous reactions are added after all the homogeneous reactions
            if any([r.reaction_type == "interface-Arrhenius",
                    r.reaction_type == "sticking-Arrhenius"]):
                continue

            if r.third_body is None:
                reaction_info.idxmap[k] = i
                reaction_info.rs.append(r)
                i += 1
    reaction_info.index.append(i)

    # everything else
    for k, r in enumerate(reaction_info.rs_unsorted):
        if r not in reaction_info.rs:
            # skip all interface-Arrhenius and sticking-Arrhenius heterogeneous reactions
            # heterogeneous reactions are added after all the homogeneous reactions
            if any([r.reaction_type == "interface-Arrhenius",
                    r.reaction_type == "sticking-Arrhenius"]):
                continue

            reaction_info.idxmap[k] = i
            reaction_info.rs.append(r)
            i += 1
    reaction_info.index.append(i)

    # surface reactions
    if not isinstance(interface, type(None)):
        for k, r in enumerate(reaction_info.rs_unsorted):
            if r not in reaction_info.rs:
                if r.reaction_type == "interface-Arrhenius":
                    reaction_info.idxmap[k] = i
                    reaction_info.rs.append(r)
                    i += 1
        reaction_info.index.append(i)

        for k, r in enumerate(reaction_info.rs_unsorted):
            if r not in reaction_info.rs:
                if r.reaction_type == "sticking-Arrhenius":
                    reaction_info.idxmap[k] = i
                    reaction_info.rs.append(r)
                    i += 1
        reaction_info.index.append(i)

    reaction_info.is_sorted = True

    return reaction_info


def rmap(fstream, mechanism, reaction_info):
    """Write reverse reaction map."""
    rmap = reaction_info.idxmap.keys()
    n_reactions = mechanism.n_reactions
    str_rmap = ",".join(str(x) for x in rmap)
    if n_reactions > 0:
        cw.writer(fstream, f"const int rmap[{n_reactions}] = {{{str_rmap}}};")


def get_rmap(fstream, mechanism):
    """Write function for reverse reaction map."""
    n_reactions = mechanism.n_reactions
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns 0-based map of reaction order"))
    cw.writer(fstream, "void GET_RMAP" + cc.sym)

    if n_reactions > 0:
        cw.writer(fstream, "(int * _rmap)")
    else:
        cw.writer(fstream, "(int * /*_rmap*/)")

    cw.writer(fstream, "{")

    if n_reactions > 0:
        cw.writer(fstream, f"for (int j=0; j<{n_reactions}; ++j)")
        cw.writer(fstream, "{")
        cw.writer(fstream, "_rmap[j] = rmap[j];")
        cw.writer(fstream, "}")

    cw.writer(fstream, "}")
