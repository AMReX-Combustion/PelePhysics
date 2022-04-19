import sys
from collections import OrderedDict

import cantera as ct


class ReactionInfo:
    """Information on reactions."""

    def __init__(self, mechanism):
        self.nReactions = mechanism.n_reactions
        self.rs = []
        self.rs_unsorted = mechanism.reactions()
        self.sorted = False
        self.index = [0]
        self.idxmap = OrderedDict()
        self.reacRemoveIDList = []  # FIXME


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

    reaction_info.sorted = True

    return reaction_info
