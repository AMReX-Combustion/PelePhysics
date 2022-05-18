"""Reaction information."""
from collections import OrderedDict

import ceptr.constants as cc
import ceptr.writer as cw


class ReactionInfo:
    """Information on reactions."""

    def __init__(self, mechanism):
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


def sort_reactions(mechanism):
    """Sort reactions."""
    reaction_info = ReactionInfo(mechanism)
    i = 0

    # troe
    for k, r in enumerate(reaction_info.rs_unsorted):
        if r not in reaction_info.rs:
            if r.reaction_type == "falloff":
                if r.rate.type == "Troe":
                    reaction_info.idxmap[k] = i
                    reaction_info.rs.append(r)
                    i += 1
    reaction_info.index.append(i)

    # sri
    for k, r in enumerate(reaction_info.rs_unsorted):
        if r not in reaction_info.rs:
            if r.reaction_type == "falloff":
                if r.rate.type == "Sri":
                    reaction_info.idxmap[k] = i
                    reaction_info.rs.append(r)
                    i += 1
    reaction_info.index.append(i)

    # lindemann
    for k, r in enumerate(reaction_info.rs_unsorted):
        if r not in reaction_info.rs:
            if r.reaction_type == "falloff":
                if r.rate.type == "Lindemann":
                    reaction_info.idxmap[k] = i
                    reaction_info.rs.append(r)
                    i += 1
    reaction_info.index.append(i)

    # three-body
    for k, r in enumerate(reaction_info.rs_unsorted):
        if r not in reaction_info.rs:
            if r.reaction_type == "three-body":
                reaction_info.idxmap[k] = i
                reaction_info.rs.append(r)
                i += 1
    reaction_info.index.append(i)

    # simplest case
    for k, r in enumerate(reaction_info.rs_unsorted):
        if r not in reaction_info.rs:
            if r.reaction_type != "three-body":
                reaction_info.idxmap[k] = i
                reaction_info.rs.append(r)
                i += 1
    reaction_info.index.append(i)

    # everything else
    for k, r in enumerate(reaction_info.rs_unsorted):
        if r not in reaction_info.rs:
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
    cw.writer(fstream, "const int rmap[%d] = {%s};" % (n_reactions, str_rmap))


def get_rmap(fstream, mechanism):
    """Write function for reverse reaction map."""
    n_reactions = mechanism.n_reactions
    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Returns 0-based map of reaction order"))
    cw.writer(fstream, "void GET_RMAP" + cc.sym + "(int * _rmap)")
    cw.writer(fstream, "{")

    cw.writer(fstream, "for (int j=0; j<%d; ++j) {" % (n_reactions))
    cw.writer(fstream, "_rmap[j] = rmap[j];")
    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
