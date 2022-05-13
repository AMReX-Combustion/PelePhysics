"""QSSA functions needed for conversion."""
import sys
from collections import OrderedDict

import numpy as np

import ceptr.utilities as cu
import ceptr.writer as cw


def set_qssa_reactions(mechanism, species_info, reaction_info):
    """Get list of reaction indices that involve QSSA species."""
    for orig_idx, idx in reaction_info.idxmap.items():
        reaction = mechanism.reaction(orig_idx)
        qssa_reaction = False
        lst_reactants = [(k, v) for k, v in reaction.reactants.items()]
        lst_products = [(k, v) for k, v in reaction.products.items()]
        agents = list(set(lst_reactants + lst_products))
        for symbol, coefficient in agents:
            if symbol in species_info.qssa_species_list:
                qssa_reaction = True
        if qssa_reaction:
            reaction_info.qssa_reactions.append(orig_idx)
            print(
                "Found a QSSA reaction: ",
                orig_idx,
                reaction.equation,
            )
    reaction_info.n_qssa_reactions = len(reaction_info.qssa_reactions)


def get_qssa_networks(mechanism, species_info, reaction_info):
    """QSSA networks.

    Get networks showing which QSSA species are involved with one
    another and which reactions each QSSA species is involved in

    """
    # Create QSSA networks for mechanism
    create_ss_net(mechanism, species_info, reaction_info)
    create_sr_net(mechanism, species_info, reaction_info)

    # Get non-zero indices of networks to be used for coupling
    # "i" is for row "j" is for column
    species_info.qssa_ss_si, species_info.qssa_ss_sj = np.nonzero(
        species_info.qssa_ssnet
    )
    species_info.sr_si, species_info.qssa_sr_rj = np.nonzero(
        species_info.qssa_srnet
    )

    print("\n\n SS network for QSSA: ")
    print(species_info.qssa_ssnet)
    print(" SR network for QSSA: ")
    print(species_info.qssa_srnet)


def create_ss_net(mechanism, species_info, reaction_info):
    """Create the species-species network."""
    for orig_idx, idx in reaction_info.idxmap.items():
        reaction = mechanism.reaction(orig_idx)

        remove_forward = cu.is_remove_forward(reaction_info, orig_idx)

        slist = []
        # get a list of species involved in the reactants and products
        for symbol, coefficient in reaction.reactants.items():
            if symbol in species_info.qssa_species_list:
                slist.append(symbol)
        for symbol, coeffecient in reaction.products.items():
            if symbol in species_info.qssa_species_list:
                slist.append(symbol)

        # if species s1 and species s2 are in the same reaction,
        # denote they are linked in the species-species network
        for s1 in slist:
            for s2 in slist:
                # we should not use the original indices, but the reordered one
                # species_info.qssa_ssnet[mechanism.qssa_species(s1).id][mechanism.qssa_species(s2).id] = 1
                species_info.qssa_ssnet[
                    species_info.ordered_idx_map[s1] - species_info.n_species
                ][
                    species_info.ordered_idx_map[s2] - species_info.n_species
                ] = 1


def create_sr_net(mechanism, species_info, reaction_info):
    """Create the species-reac network."""
    # for each reaction in the mechanism
    for orig_idx, idx in reaction_info.idxmap.items():
        reaction = mechanism.reaction(orig_idx)
        reactant_list = []
        product_list = []

        remove_forward = cu.is_remove_forward(reaction_info, orig_idx)

        # get a list of species involved in the reactants and products

        for symbol, coefficient in reaction.reactants.items():
            if symbol in species_info.qssa_species_list:
                reactant_list.append(symbol)
        for symbol, coeffecient in reaction.products.items():
            if symbol in species_info.qssa_species_list:
                product_list.append(symbol)

        # if qss species s is in reaction number i,
        # denote they are linked in the Species-Reaction network
        # with negative if s is a reactant(consumed) and positive if s is a product(produced)
        for s in reactant_list:
            species_info.qssa_srnet[
                species_info.ordered_idx_map[s] - species_info.n_species
            ][orig_idx] = -1
        for s in product_list:
            species_info.qssa_srnet[
                species_info.ordered_idx_map[s] - species_info.n_species
            ][orig_idx] = 1


def qssa_validation(mechanism, species_info, reaction_info):
    """Check that the QSSA given are actually "valid" options.

    Exit if species is not valid.
    """
    # Check that QSSA species are all species used in the given mechanism
    for s in species_info.qssa_species_list:
        if s not in species_info.all_species_list:
            text = "species " + s + " is not in the mechanism"
            sys.exit(text)

    # Check that QSSA species are consumed/produced at least once to ensure theoretically valid QSSA option
    # (There is more to it than that, but this is a quick catch based on that aspect)
    for i, symbol in enumerate(species_info.qssa_species_list):
        consumed = 0
        produced = 0
        for orig_idx in species_info.qssa_sr_rj[species_info.sr_si == i]:
            reaction = mechanism.reaction(orig_idx)
            remove_forward = cu.is_remove_forward(reaction_info, orig_idx)
            if any(
                reactant == symbol
                for reactant, _ in reaction.reactants.items()
            ):
                if remove_forward:
                    consumed += 0
                else:
                    consumed += 1
                if reaction.reversible:
                    produced += 1
            if any(
                product == symbol for product, _ in reaction.products.items()
            ):
                if remove_forward:
                    produced += 0
                else:
                    produced += 1
                if reaction.reversible:
                    consumed += 1

        if consumed == 0 or produced == 0:
            text = (
                "Uh Oh! QSSA species "
                + symbol
                + " does not have a balanced consumption/production relationship in mechanism => bad QSSA choice"
            )
            sys.exit(text)


def qssa_coupling(mechanism, species_info, reaction_info):
    """Determine from qssa_SSnet which QSSA species depend on each other specifically."""
    is_coupling = False
    coupling_reactions = ""
    list_coupling_reactions = []

    species_info.qssa_scnet = species_info.qssa_ssnet
    for i in range(species_info.n_qssa_species):
        for j in species_info.qssa_ss_sj[species_info.qssa_ss_si == i]:
            if j != i:
                count = 0
                for r in species_info.qssa_sr_rj[species_info.sr_si == j]:
                    reaction = mechanism.reaction(r)

                    remove_forward = cu.is_remove_forward(reaction_info, r)

                    # put forth any pathological case
                    if any(
                        reactant == species_info.qssa_species_list[j]
                        for reactant, _ in reaction.reactants.items()
                    ):
                        if any(
                            product == species_info.qssa_species_list[j]
                            for product, _ in reaction.products.items()
                        ):
                            sys.exit(
                                "Species "
                                + species_info.qssa_species_list[j]
                                + " appears as both prod and reacts. Check reaction "
                                + reaction.equation()
                            )

                    # we know j is in reaction r. Options are
                    # IF r is reversible
                    # j + i <-> prods OR reacts <-> j + i NOT ALLOWED
                    # all other combinations are fine.
                    # IF r is not reversible
                    # j + i -> prods NOT ALLOWED
                    # all other combinations are fine
                    # note that it is assumed no coupling bet same QSSA -- this is if i == j

                    if reaction.reversible:
                        # QSSA spec j is a reactant
                        if any(
                            reactant == species_info.qssa_species_list[j]
                            for reactant, _ in reaction.reactants.items()
                        ):
                            # Check if QSSA species i is a reactant too
                            if (not remove_forward) and any(
                                reactant == species_info.qssa_species_list[i]
                                for reactant, _ in reaction.reactants.items()
                            ):
                                is_coupling = True
                                if reaction.id not in list_coupling_reactions:
                                    coupling_reactions += (
                                        "R"
                                        + str(reaction.orig_id)
                                        + " "
                                        + reaction.equation()
                                        + "\n"
                                    )
                                    list_coupling_reactions.append(reaction.id)
                                print(
                                    "Quadratic coupling between "
                                    + species_info.qssa_species_list[j]
                                    + " and "
                                    + species_info.qssa_species_list[i]
                                    + " in reaction "
                                    + reaction.equation()
                                    + " not allowed !!!"
                                )
                            # if QSSA specices j is a reactant and QSSA species i is a product,
                            # because react is two way then j depend on i and vice-versa
                            elif any(
                                product == species_info.qssa_species_list[i]
                                for product, _ in reaction.products.items()
                            ):
                                count += 1
                        # if QSSA species j is not a reactant, then it must be a product.
                        else:
                            # Check if QSSA species i is also a product
                            if any(
                                product == species_info.qssa_species_list[i]
                                for product, _ in reaction.products.items()
                            ):
                                is_coupling = True
                                if reaction.id not in list_coupling_reactions:
                                    coupling_reactions += (
                                        "R"
                                        + str(reaction.orig_id)
                                        + " "
                                        + reaction.equation()
                                        + "\n"
                                    )
                                    list_coupling_reactions.append(reaction.id)
                                print(
                                    "Quadratic coupling between "
                                    + species_info.qssa_species_list[j]
                                    + " and "
                                    + species_info.qssa_species_list[i]
                                    + " in reaction "
                                    + reaction.equation()
                                    + " not allowed !!!"
                                )
                            # if QSSA specices j is a product and QSSA species i is a reactant
                            # because react is two way then j depend on i and vice-versa
                            elif any(
                                reactant == species_info.qssa_species_list[i]
                                for reactant, _ in reaction.reactants.items()
                            ):
                                count += 1
                    else:
                        # QSSA spec j is a reactant
                        if any(
                            reactant == species_info.qssa_species_list[j]
                            for reactant, _ in reaction.reactants.items()
                        ):
                            # Check if QSSA species i is a reactant too
                            if any(
                                reactant == species_info.qssa_species_list[i]
                                for reactant, _ in reaction.reactants.items()
                            ):
                                is_coupling = True
                                if reaction.id not in list_coupling_reactions:
                                    coupling_reactions += (
                                        "R"
                                        + str(reaction.orig_id)
                                        + " "
                                        + reaction.equation()
                                        + "\n"
                                    )
                                    list_coupling_reactions.append(reaction.id)
                                print(
                                    "Quadratic coupling between "
                                    + species_info.qssa_species_list[j]
                                    + " and "
                                    + species_info.qssa_species_list[i]
                                    + " in reaction "
                                    + reaction.equation()
                                    + " not allowed !!!"
                                )
                            # if QSSA specices j is a reactant and QSSA species i is a product
                            elif any(
                                product == species_info.qssa_species_list[i]
                                for product, _ in reaction.products.items()
                            ):
                                count += 1

                if count == 0:
                    # i depends on j
                    species_info.qssa_scnet[i, j] = 0
            else:
                species_info.qssa_scnet[i, j] = 0
                for r in species_info.qssa_sr_rj[species_info.sr_si == j]:
                    reaction = mechanism.reaction(r)
                    remove_forward = cu.is_remove_forward(reaction_info, r)
                    species_appearances = 0

                    # put forth any pathological case
                    if any(
                        reactant == species_info.qssa_species_list[j]
                        for reactant, _ in reaction.reactants.items()
                    ):
                        if any(
                            product == species_info.qssa_species_list[j]
                            for product, _ in reaction.products.items()
                        ):
                            sys.exit(
                                "Species "
                                + species_info.qssa_species_list[j]
                                + " appears as both prod and reacts. Check reaction "
                                + reaction.equation()
                            )

                    if reaction.reversible:
                        # QSSA j is a reactant
                        if (not remove_forward) and any(
                            reactant == species_info.qssa_species_list[j]
                            for reactant, _ in reaction.reactants.items()
                        ):
                            for spec, coeff in reaction.reactants.items():
                                if spec == species_info.qssa_species_list[j]:
                                    species_appearances += 1
                                    if (coeff > 1.0) or (
                                        species_appearances > 1
                                    ):
                                        is_coupling = True
                                        if (
                                            reaction.id
                                            not in list_coupling_reactions
                                        ):
                                            coupling_reactions += (
                                                "R"
                                                + str(reaction.orig_id)
                                                + " "
                                                + reaction.equation()
                                                + "\n"
                                            )
                                            list_coupling_reactions.append(
                                                reaction.id
                                            )
                                        print(
                                            "Quadratic coupling of "
                                            + species_info.qssa_species_list[j]
                                            + " with itself in reaction "
                                            + reaction.equation()
                                            + " not allowed !!!"
                                        )
                        # if QSSA species j is not a reactant, then it must be a product.
                        else:
                            for spec, coeff in reaction.products.items():
                                if spec == species_info.qssa_species_list[j]:
                                    species_appearances += 1
                                    if (coeff > 1.0) or (
                                        species_appearances > 1
                                    ):
                                        is_coupling = True
                                        if (
                                            reaction.id
                                            not in list_coupling_reactions
                                        ):
                                            coupling_reactions += (
                                                "R"
                                                + str(reaction.orig_id)
                                                + " "
                                                + reaction.equation()
                                                + "\n"
                                            )
                                            list_coupling_reactions.append(
                                                reaction.id
                                            )
                                        print(
                                            "Quadratic coupling of "
                                            + species_info.qssa_species_list[j]
                                            + " with itself in reaction "
                                            + reaction.equation()
                                            + " not allowed !!!"
                                        )
                    else:
                        # QSSA spec j is a reactant
                        if any(
                            reactant == species_info.qssa_species_list[j]
                            for reactant, _ in reaction.reactants.items()
                        ):
                            for spec, coeff in reaction.reactants.items():
                                if spec == species_info.qssa_species_list[j]:
                                    species_appearances += 1
                                    if (coeff > 1.0) or (
                                        species_appearances > 1
                                    ):
                                        is_coupling = True
                                        if (
                                            reaction.id
                                            not in list_coupling_reactions
                                        ):
                                            coupling_reactions += (
                                                "R"
                                                + str(reaction.orig_id)
                                                + " "
                                                + reaction.equation()
                                                + "\n"
                                            )
                                            list_coupling_reactions.append(
                                                reaction.id
                                            )
                                        print(
                                            "Quadratic coupling of "
                                            + species_info.qssa_species_list[j]
                                            + " with itself in reaction "
                                            + reaction.equation()
                                            + " not allowed !!!"
                                        )

    species_info.qssa_sc_si, species_info.qssa_sc_sj = np.nonzero(
        species_info.qssa_scnet
    )
    print("\n\n SC network for QSSA: ")
    print(species_info.qssa_scnet)
    if is_coupling:
        sys.exit(
            "There is some quadratic coupling in mechanism. Here is the list of reactions to check: \n"
            + coupling_reactions
        )


def set_qssa_needs(mechanism, species_info, reaction_info):
    """QSSA needs."""
    self.needs = OrderedDict()
    self.needs_count = OrderedDict()

    self.needs_running = OrderedDict()
    self.needs_count_running = OrderedDict()

    for i in range(species_info.n_qssa_species):
        needs_species = []
        count = 0
        for j in species_info.qssa_sc_sj[species_info.qssa_sc_si == i]:
            if j != i:
                needs_species.append(species_info.qssa_species_list[j])
                count += 1
        self.needs[species_info.qssa_species_list[i]] = needs_species
        self.needs_count[species_info.qssa_species_list[i]] = count

    self.needs_running = self.needs.copy()
    self.needs_count_running = self.needs_count.copy()

    print("Needs report (one per QSSA spec): ")
    print(self.needs)


def set_qssa_isneeded(mechanism, species_info, reaction_info):
    """QSSA is needed."""
    self.is_needed = OrderedDict()
    self.is_needed_count = OrderedDict()

    self.is_needed_running = OrderedDict()
    self.is_needed_count_running = OrderedDict()

    for i in range(species_info.n_qssa_species):
        is_needed_species = []
        count = 0
        for j in species_info.qssa_sc_si[species_info.qssa_sc_sj == i]:
            if j != i:
                is_needed_species.append(species_info.qssa_species_list[j])
                count += 1
        self.is_needed[species_info.qssa_species_list[i]] = is_needed_species
        self.is_needed_count[species_info.qssa_species_list[i]] = count

    self.is_needed_running = self.is_needed.copy()
    self.is_needed_count_running = self.is_needed_count.copy()

    print("Is needed report (one per QSSA spec): ")
    print(self.is_needed)


def get_qssa_groups(mechanism, species_info, reaction_info):
    """Get strongly connected components (groups).

    this includes two-way dependencies: (s1 needs s2) and (s2 needs s1) => group
    and cyclical dependencies: (s1 needs s2) and (s2 needs s3) and (s3 needs s1) => group

    This function along with find_closed_cycle is an implmentation of:

    Tarjan's Strongly Connected Components Algorithm (1972)

    where node "index" is stored implictly via the location of the
    node in the lowest_link OrderedDict and boolean "onStack"
    information is stored implicitly by checking for node existence in
    potential_group

    """
    self.group = OrderedDict()

    self.discovery_order = 0

    self.potential_group = []
    self.lowest_link = OrderedDict()

    print("\n\nDetermining groups of coupled species now...")
    print("---------------------------------")
    self.all_groups = OrderedDict()

    # Loop through species to tackle the needs group
    for member in self.needs_running.keys():
        # Only need to check things that have not already been searched; i.e. they have no id/don't exist in the lowest link value list yet
        if member not in self.lowest_link.keys():

            print("- dealing with group: ", member)
            find_closed_cycle(mechanism, member)

    print("** Groups of coupled species are: ", self.all_groups)

    group_count = 0
    # Rename and only store strongly connected components involving 2 or more species as groups
    for group in self.all_groups:
        if len(self.all_groups[group]) > 1:
            self.group["group_" + str(group_count)] = self.all_groups[group]
            group_count += 1
    print()
    print("** Final clean self groups are: ", self.group)

    update_group_needs(mechanism, species_info, reaction_info)
    update_group_dependencies(mechanism, species_info, reaction_info)


def find_closed_cycle(mechanism, species):
    """Find closed cycle."""

    # We start the recursion on node "species".
    # This species is considered the "parent" node.
    print(
        "      Searching for potential closed cycle involving parent: ",
        species,
    )

    # We only enter the recursion if the species has not already been discovered
    # so we add this species to the potential group list and give it a lowest link value
    # The species location in the lowest link value dictionary denotes its discovery order, or id
    self.potential_group.append(species)
    self.lowest_link[species] = self.discovery_order
    self.discovery_order += 1
    parent = species

    print(
        "      Upon initialization, discovery order is: ",
        list(self.lowest_link.keys()).index(species),
        " and lowest link value is: ",
        self.discovery_order - 1,
    )
    print()

    # Loop through the needs of the parent
    # Each need is denoted as a "child" of the parent
    for need in self.needs_running[parent]:
        child = need
        print("       x Start level of needs loop ")
        print("       x Child is: ", child)

        # If the child has not yet been discovered, we recurse so that the child becomes the parent and the search continues as described above
        if child not in self.lowest_link.keys():

            print("         xx Child has never been visited at all...")
            print("         xx Initiate recursion to assign lowlink value ")

            find_closed_cycle(mechanism, child)

            print(
                "         xx We've finished a recursion! The child that was passed in was: ",
                child,
                " with parent ",
                parent,
            )
            # print("         The discovery order of the parent is: ", list(self.lowest_link.keys()).index(parent))
            print(
                "         xx The lowest link value of the parent is: ",
                self.lowest_link[parent],
            )
            # print("         The discovery order of this child is: ", list(self.lowest_link.keys()).index(child))
            print(
                "         xx The lowest link value of this child is: ",
                self.lowest_link[child],
            )
            print()

            # If the lowest link connection of the child is lower than that of the parent's, then this lowest link value becomes the parent's new lowest link value
            # This update comes into effect when the child at the end of the recursion
            # matches the original parent of the search or has a lower link connection from an earlier search than the parent currently does
            if child in self.potential_group:
                self.lowest_link[parent] = min(
                    self.lowest_link[parent], self.lowest_link[child]
                )

            print(
                "         ** Therefore, the lowest link value of the parent then becomes: ",
                self.lowest_link[parent],
            )
            print(
                "         ** The lowest link connections are now: ",
                self.lowest_link,
            )
            print()
            print()

        # If the child is already listed in the potential_group, then it has been discovered during this search and is already a part of the simply connected component
        # Note: If the child already has a lowest link value but is NOT in the potential_group, that means it is a part of a previously found group
        elif child in self.potential_group:

            print("         xx Child has already been visited this recursion ")
            # print("         The discovery order of the parent is: ", list(self.lowest_link.keys()).index(parent))
            print(
                "         xx The lowest link value of the parent is: ",
                self.lowest_link[parent],
            )
            print(
                "         xx The discovery order of this child is: ",
                list(self.lowest_link.keys()).index(child),
            )
            # print("         The lowest link value of this child is: ", self.lowest_link[child])

            # Since the child has been discovered already during this search, that means it is in the group but still in the recursion process,
            # Update the parent's lowest link value with this childs discovery order
            self.lowest_link[parent] = min(
                self.lowest_link[parent],
                list(self.lowest_link.keys()).index(child),
            )

            print()
            print(
                "         **** Therefore, the lowest link value of the parent then becomes: ",
                self.lowest_link[parent],
            )
            print(
                "         **** The lowest link connections are now: ",
                self.lowest_link,
            )
            print()
            print()

    # If, after searching all children and updating lowest link connections, you reach a parent whose lowest link still matches its original value (or starting position)
    # Then you have reached the root of a simply connected component (aka you've found a group)
    if list(self.lowest_link.keys()).index(parent) == self.lowest_link[parent]:
        hold = []
        while True:
            # remove group species from the running potential group list until you reach the parent species
            # add these to the official group, leaving other potential group components in the potential_group list for continued recursion completion
            node = self.potential_group.pop()
            hold.append(node)
            if node == parent:
                break
        self.all_groups[parent] = hold
        print("         Group is: ", self.all_groups[parent])
        print()
        print()


def update_group_needs(mechanism, species_info, reaction_info):
    """Update group member needs with group names.

    group member needs species -> group needs species
    group member is needed by species -> group is needed by species
    """
    print("\n\nUpdating group needs...")
    print("---------------------------------")

    for group_key in list(self.group.keys()):
        print("-Dealing with group " + group_key)

        update_needs = []
        update_is_needed = []
        update_needs_count = 0
        update_needed_count = 0

        group_needs = {}
        group_needs_count = {}
        group_is_needed = {}
        group_is_needed_count = {}

        other_groups = list(self.group.keys())
        other_groups.remove(group_key)
        print("  (other groups are: ", other_groups, ")")

        # for each species in the current group
        for spec in self.group[group_key]:
            print("... for group member: " + spec)
            # look at any additional needs that are not already accounted for with the group
            for need in list(
                set(self.needs_running[spec]) - set(self.group[group_key])
            ):
                print("        An additional not-in-group need is " + need)
                not_in_group = True
                # check the other groups to see if the need can be found in one of them
                for other_group in other_groups:
                    # if the other group is not already accounted for
                    # and it contains the spec need we're looking for,
                    # update the group needs with that group that contains the spec need
                    if other_group not in update_needs and any(
                        member == need for member in self.group[other_group]
                    ):
                        print(
                            "        it is found in a different group. Adding it."
                        )
                        not_in_group = False
                        update_needs.append(other_group)
                        update_needs_count += 1
                    elif other_group in update_needs and any(
                        member == need for member in self.group[other_group]
                    ):
                        print(
                            "        it is found in a group that was already put in the list due to the fact that another species in the group is needed by the current species."
                        )
                        not_in_group = False
                # alternatively, if this is just a solo need that's not in another group,
                # update the group needs with just that need.
                if not_in_group and need not in update_needs:
                    print(
                        "        this need was not found in a group ! Adding the spec directly"
                    )
                    update_needs.append(need)
                    update_needs_count += 1
            # look at any additional species (outside of the group) that depend on the current group member
            for needed in list(
                set(self.is_needed_running[spec]) - set(self.group[group_key])
            ):
                print(
                    "        An additional not-in-group is-needed is " + needed
                )
                not_in_group = True
                # for the other groups
                for other_group in other_groups:
                    # if the other group hasn't alredy been accounted for and the species is in that group, then that other group depends on a species in the current group
                    if other_group not in update_is_needed and any(
                        member == needed for member in self.group[other_group]
                    ):
                        print(
                            "        it is found in a different group. Adding it."
                        )
                        not_in_group = False
                        update_is_needed.append(other_group)
                        update_needed_count += 1
                    elif other_group in update_is_needed and any(
                        member == needed for member in self.group[other_group]
                    ):
                        print(
                            "        it is foud in a group that was already put in the list due to the fact that another species in the group is needed by the current species."
                        )
                        not_in_group = False
                # if the species is not in another group, then that lone species just depends on the current group.
                if not_in_group and needed not in update_is_needed:
                    print(
                        "        this is-needed was not found in a group ! Adding the spec directly"
                    )
                    update_is_needed.append(needed)
                    update_needed_count += 1

            # del self.needs_running[spec]
            # del self.needs_count_running[spec]
            # del self.is_needed_running[spec]

        group_needs[group_key] = update_needs
        group_needs_count[group_key] = update_needs_count
        group_is_needed[group_key] = update_is_needed
        group_is_needed_count[group_key] = update_needed_count

        self.needs_running.update(group_needs)
        self.needs_count_running.update(group_needs_count)
        self.is_needed_running.update(group_is_needed)
        self.is_needed_count_running.update(group_is_needed_count)

        print("So, ", group_key, " needs ", update_needs)
        print("So, ", group_key, " is-needed is ", update_is_needed)

    for group in list(self.group.keys()):
        for spec in self.group[group]:
            if spec in self.needs_running:
                del self.needs_running[spec]
                del self.needs_count_running[spec]
                del self.is_needed_running[spec]

    print()
    print("** This is the final needs running and is_needed running: ")
    print(self.needs_running)
    print(self.is_needed_running)


def update_group_dependencies(mechanism, species_info, reaction_info):
    """Update solo species dependendent on group members with group names.

    species needs member -> species needs group
    species is needed by group member -> species is needed by group
    """
    print("\n\nUpdating group dependencies...")
    print("---------------------------------")

    solo_needs = self.needs_running.copy()
    solo_needs_count = self.needs_count_running.copy()
    solo_is_needed = self.is_needed_running.copy()
    solo_is_needed_count = self.is_needed_count_running.copy()

    # remove the groups because we're just dealing with things that aren't in groups now
    for group in list(self.group.keys()):
        del solo_needs[group]
        del solo_needs_count[group]
        del solo_is_needed[group]
        del solo_is_needed_count[group]

    for solo in list(solo_needs.keys()):
        print("-Dealing with solo species " + solo)
        update_needs = []
        update_is_needed = []
        update_needs_count = 0
        update_needed_count = 0
        for need in solo_needs[solo]:
            print("... who needs: " + need)
            not_in_group = True
            for group in list(self.group.keys()):
                if group not in update_needs and any(
                    member == need for member in self.group[group]
                ):
                    print("        this species is in group: ", group)
                    not_in_group = False
                    update_needs.append(group)
                    update_needs_count += 1
                elif group in update_needs and any(
                    member == need for member in self.group[group]
                ):
                    print(
                        "        this group was already put in the list due to the fact that another species in the group is needed by the current species."
                    )
                    not_in_group = False
            if not_in_group and need not in update_needs:
                print(
                    "        this need was not found in a group ! Adding the spec directly"
                )
                update_needs.append(need)
                update_needs_count += 1

        for needed in solo_is_needed[solo]:
            print("... who is-needed needs are: " + needed)
            not_in_group = True
            for group in list(self.group.keys()):
                if group not in update_is_needed and any(
                    member == needed for member in self.group[group]
                ):
                    print("        this species is in group: ", group)
                    not_in_group = False
                    update_is_needed.append(group)
                    update_needed_count += 1
                if group in update_is_needed and any(
                    member == needed for member in self.group[group]
                ):
                    print(
                        "        this group was already put in the list due to the fact that another species in the group is needed by the current species."
                    )
                    not_in_group = False
            if not_in_group and needed not in update_is_needed:
                print(
                    "        this is-needed need was not found in a group ! Adding the spec directly"
                )
                update_is_needed.append(needed)
                update_needed_count += 1

        solo_needs[solo] = update_needs
        solo_needs_count[solo] = update_needs_count
        solo_is_needed[solo] = update_is_needed
        solo_is_needed_count[solo] = update_needed_count

    self.needs_running.update(solo_needs)
    self.needs_count_running.update(solo_needs_count)
    self.is_needed_running.update(solo_is_needed)
    self.is_needed_count_running.update(solo_is_needed_count)

    print()
    print("** This is the final needs running and is_needed running: ")
    print(self.needs_running)
    print(self.is_needed_running)


def sort_qssa_computation(mechanism, species_info, reaction_info):
    """Sort order that QSSA species need to be computed based on dependencies."""
    self.decouple_index = OrderedDict()
    self.decouple_count = 0

    # look at how many dependencies each component has
    needs_count_regress = self.needs_count_running.copy()

    # There should always be a component present that loses
    # all dependencies as you update the computation
    while 0 in list(needs_count_regress.values()):
        needs_count_base = needs_count_regress.copy()
        # for each component (species, group, sup group, etc.) that needs things...
        for member in needs_count_base:
            # if that component doesn't need anything
            if needs_count_base[member] == 0:
                print("-Dealing with member ", member)
                # solve that component now
                self.decouple_index[self.decouple_count] = member
                # then delete it out of the updating needs list
                del needs_count_regress[member]
                # for anything needed by that component
                for needed in self.is_needed_running[member]:
                    # decrease that thing's dependency since it has now been taken care of
                    needs_count_regress[needed] -= 1
                self.decouple_count += 1

    # If your decouple count doesn't match the number of components with needs,
    # then the system is more complicated than what these functions can handle currently
    if len(self.decouple_index) != len(self.needs_running):
        print("WARNING: Some components may not have been taken into account")
    print()
    print(
        "** order of execution for qss concentration calculations: ",
        self.decouple_index,
    )


def sort_qssa_solution_elements(mechanism, species_info, reaction_info):
    """Components needed to set up QSSA algebraic expressions from AX = B.

    where A contains coefficients from qf's and qr's, X contains QSSA species concentrations,
    and B contains qf's and qr's
    Info stored as: RHS vector (non-QSSA and QSSA qf's and qr's),
    coefficient of species (diagonal elements of A), coefficient of group mates (coupled off-diagonal elements of A)
    """
    self.qssa_rhs = OrderedDict()
    self.qssa_coeff = OrderedDict()
    # self.qssa_groupSp   = OrderedDict()
    self.qssa_qssa_coeff = OrderedDict()

    # Need to get qfqr_coeff reaction map
    ispecial = reaction_info.index[5:7]
    nspecial_qss = 0
    ispecial_qss = [0, 0]
    special_first = True

    # Find out bounds for special reacs in smaller qssReactions list
    for reac_id in self.qssReactions:
        if reac_id >= ispecial[0] and reac_id < ispecial[1]:
            nspecial_qss += 1
            if special_first:
                ispecial_qss[0] = self.qssReactions.index(reac_id)
                special_first = False
            ispecial_qss[1] = self.qssReactions.index(reac_id) + 1

    # remove special reacs for some reason ?
    self.qfqr_co_idx_map = self.qssReactions
    if (ispecial_qss[1] - ispecial_qss[0]) > 0:
        for index in range(ispecial_qss[0], ispecial_qss[1]):
            del self.qfqr_co_idx_map[index]

    for i in range(species_info.n_qssa_species):
        symbol = species_info.qssa_species_list[i]
        print()
        print("-<>-Dealing with QSSA species ", i, symbol)
        print("__________________________________________")
        coupled = []
        reactants = []
        products = []
        rhs_hold = []
        coeff_hold = []
        groupCoeff_hold = defaultdict(list)

        for r in species_info.qssa_sr_rj[species_info.sr_si == i]:
            reaction = mechanism.reaction(id=r)
            remove_forward = cu.is_remove_forward(reaction_info, orig_idx)
            # check this mess of reactions
            if (reaction.id - 1) != r:
                print("\n\nCheck this!!!\n")
                sys.exit(1)
            print("... who is involved in reac ", r, reaction.equation())
            print(
                "... reaction ",
                reaction.id,
                "is QSSA reaction number ",
                self.qfqr_co_idx_map.index(reaction.id - 1),
            )

            direction = species_info.qssa_srnet[i][r]

            # Check if reaction contains other QSSA species
            coupled = [
                species
                for species in list(
                    set(species_info.sr_si[species_info.qssa_sr_rj == r])
                )
            ]
            coupled_qss = [species_info.qssa_species_list[j] for j in coupled]

            for species in coupled_qss:
                if species != symbol:
                    other_qss = species

            if len(coupled) >= 2:
                # Check if all QSSA are only reactants
                all_qssa_reactants = True
                if any(
                    product == other_qss
                    for product, _ in reaction.products.items()
                ):
                    all_qssa_reactants = False
                if any(
                    product == symbol
                    for product, _ in reaction.products.items()
                ):
                    all_qssa_reactants = False

            if len(coupled) < 2:
                print("        this reaction only involves that QSSA ")
                # if QSSA species is a reactant
                if direction == -1:
                    print(
                        "        species ",
                        symbol,
                        " in reaction ",
                        r,
                        " is a reactant",
                    )
                    species_appearances = 0
                    for reactant in reaction.reactants:
                        spec, coeff = reactant
                        if spec == symbol:
                            species_appearances += 1

                    coeff_hold.append(
                        "-qf_co[" + str(self.qfqr_co_idx_map.index(r)) + "]"
                    )
                    if reaction.reversible:
                        rhs_hold.append(
                            "+"
                            + str(float(species_appearances))
                            + "*qr_co["
                            + str(self.qfqr_co_idx_map.index(r))
                            + "]"
                        )
                # if QSSA species is a product
                elif direction == 1:
                    print(
                        "        species ",
                        symbol,
                        " in reaction ",
                        r,
                        " is a product",
                    )
                    species_appearances = 0
                    for product in reaction.products:
                        spec, coeff = product
                        if spec == symbol:
                            species_appearances += 1

                    rhs_hold.append(
                        "+"
                        + str(float(species_appearances))
                        + "*qf_co["
                        + str(self.qfqr_co_idx_map.index(r))
                        + "]"
                    )
                    if reaction.reversible:
                        coeff_hold.append(
                            "-qr_co["
                            + str(self.qfqr_co_idx_map.index(r))
                            + "]"
                        )

            elif len(coupled) >= 2 and all_qssa_reactants and remove_forward:
                # coupling is actually not a coupling because forward reaction is removed
                if reaction.reversible:
                    # Should always be true
                    # Check how many times species appear in reactants
                    species_appearances = 0
                    for reactant in reaction.reactants:
                        spec, coeff = reactant
                        if spec == symbol:
                            species_appearances += 1

                    rhs_hold.append(
                        "+"
                        + str(float(species_appearances))
                        + "*qr_co["
                        + str(self.qfqr_co_idx_map.index(r))
                        + "]"
                    )

            else:
                # note in this case there can only be 2 QSSA in one reac
                coupled_qss = [
                    species_info.qssa_species_list[j] for j in coupled
                ]
                print(
                    "        this reaction couples the following QSSA: ",
                    coupled_qss,
                )

                # assumes only 2 QSSA can appear in a reac now
                # other_qssa_list = []
                # for species in coupled_qss:
                #    other_qssa_list.append(species)
                # other_qssa_list.remove(symbol)
                for species in coupled_qss:
                    if species != symbol:
                        other_qss = species
                # for group in self.group:
                #    if set(coupled_qss).issubset(set(self.group[group])):
                #        "        (they are both on the same group)"
                #        group_flag = True

                # THIS is the right groupCoeff list
                # if group_flag:

                # if QSSA species is a reactant (other QSSA must be a product to be coupled to be coupled here, or quadratic coupling would have been triggered earlier)
                if direction == -1:

                    all_qssa_reactants = True
                    if any(
                        product == other_qss
                        for product, _ in reaction.products.items()
                    ):
                        all_qssa_reactants = False
                    if any(
                        product == symbol
                        for product, _ in reaction.products.items()
                    ):
                        all_qssa_reactants = False

                    remove_forward = cu.is_remove_forward(
                        reaction_info, orig_idx
                    )
                    if all_qssa_reactants and remove_forward:
                        print(
                            "        species ",
                            symbol,
                            " and species ",
                            other_qss,
                            " in purely irreversible reaction ",
                            r,
                            " are both reactants",
                        )
                        print(
                            "        this reaction does not contribute to any QSSA coefficients and is thus ignored"
                        )

                    else:

                        print(
                            "        species ",
                            symbol,
                            " in reaction ",
                            r,
                            " is a reactant",
                        )

                        species_appearances = 0
                        for reactant in reaction.reactants:
                            spec, coeff = reactant
                            if spec == symbol:
                                species_appearances += 1

                        coeff_hold.append(
                            "-qf_co["
                            + str(self.qfqr_co_idx_map.index(r))
                            + "]"
                        )
                        if reaction.reversible:
                            groupCoeff_hold[other_qss].append(
                                "+"
                                + str(float(species_appearances))
                                + "*qr_co["
                                + str(self.qfqr_co_idx_map.index(r))
                                + "]"
                            )
                # if QSSA species is a product AND other QSSA species is a reactant (not guaranteed; must check that QSSA are on opposite sides of equation)
                elif direction == 1 and any(
                    reactant == other_qss
                    for reactant, _ in reaction.reactants.items()
                ):
                    print(
                        "        species ",
                        symbol,
                        " in reaction ",
                        r,
                        " is a product",
                    )
                    print("        other qss species is ", other_qss)

                    species_appearances = 0
                    for product in reaction.products:
                        spec, coeff = product
                        if spec == symbol:
                            species_appearances += 1

                    groupCoeff_hold[other_qss].append(
                        "+"
                        + str(float(species_appearances))
                        + "*qf_co["
                        + str(self.qfqr_co_idx_map.index(r))
                        + "]"
                    )
                    if reaction.reversible:
                        coeff_hold.append(
                            "-qr_co["
                            + str(self.qfqr_co_idx_map.index(r))
                            + "]"
                        )
                # last option is that BOTH QSSA are products, but the reaction is only one way, so it doesn't matter. This is ignored in the quadratic coupling check as
                # the reverse rate would be zero and thus would not affect anything anyway.
                else:
                    print(
                        "        species ",
                        symbol,
                        " and species ",
                        other_qss,
                        " in irreversible reaction ",
                        r,
                        " are both products",
                    )
                    print(
                        "        this reaction does not contribute to any QSSA coefficients and is thus ignored"
                    )

            print()
            print(
                "After dealing with QSSA species "
                + symbol
                + " in  reaction "
                + str(r)
                + " we have the following: "
            )
            print("rhs_hold is ", rhs_hold)
            print("coeff_hold is ", coeff_hold)
            print("groupCoeff_hold is ", groupCoeff_hold)
            print()

        self.qssa_rhs[symbol] = " ".join(rhs_hold)
        self.qssa_coeff[symbol] = " ".join(coeff_hold)
        self.qssa_qssa_coeff[symbol] = OrderedDict()
        for j in range(species_info.n_qssa_species):
            if j != i:
                other_qss = species_info.qssa_species_list[j]
                if other_qss in groupCoeff_hold:
                    self.qssa_qssa_coeff[symbol][other_qss] = " ".join(
                        groupCoeff_hold[other_qss]
                    )
                else:
                    self.qssa_qssa_coeff[symbol][other_qss] = "0.0"

        # for group in self.group.keys():
        #    if any(component == symbol for component in self.group[group]):
        #        self.qssa_groupSp[symbol] = groupCoeff_hold

        print("Here is everything: ")
        print()
        print("rhs: ", self.qssa_rhs[symbol])
        print("self: ", self.qssa_coeff[symbol])
        print("coupling: ", self.qssa_qssa_coeff[symbol])
        print()

    # for species in self.qssa_groupSp.keys():
    #    for coeff in self.qssa_groupSp[species].keys():
    #        self.qssa_groupSp[species][coeff] =" ".join(self.qssa_groupSp[species][coeff])

    # for symbol in self.group.keys():
    #    for s1 in self.group[symbol]:
    #        for s2 in self.group[symbol]:
    #            if s2 != s1 and not self.qssa_groupSp[s1][s2]:
    #                self.qssa_groupSp[s1][s2] = str(0.0)

    print()
    print()
    print("Final lists: ")
    print("-------------")
    print("rhs: ", self.qssa_rhs)
    print("self: ", self.qssa_coeff)
    print("coupling: ", self.qssa_qssa_coeff)
    print()
    print()


def qssa_component_functions(fstream, mechanism, species_info, reaction_info):
    """QSSA component functions."""
    itroe = reaction_info.index[0:2]
    isri = reaction_info.index[1:3]
    ilindemann = reaction_info.index[2:4]
    i3body = reaction_info.index[3:5]
    isimple = reaction_info.index[4:6]
    ispecial = reaction_info.index[5:7]

    print("troe index range is: ", itroe)
    print("sri index range is: ", isri)
    print("lindemann index range is: ", ilindemann)
    print("3body index range is: ", i3body)
    print("simple index range is: ", isimple)
    print("special index range is: ", ispecial)

    ntroe_qss = 0
    nsri_qss = 0
    nlindemann_qss = 0
    n3body_qss = 0
    nsimple_qss = 0
    nspecial_qss = 0

    itroe_qss = [0, 0]
    isri_qss = [0, 0]
    ilindemann_qss = [0, 0]
    i3body_qss = [0, 0]
    isimple_qss = [0, 0]
    ispecial_qss = [0, 0]

    troe_first = True
    sri_first = True
    lindemann_first = True
    threebody_first = True
    simple_first = True
    special_first = True

    for reac_id in self.qssReactions:
        if reac_id >= itroe[0] and reac_id < itroe[1]:
            print(
                "reaction ",
                reac_id,
                mechanism.reaction(id=reac_id).equation(),
                " goes in troe",
            )
            ntroe_qss += 1
            if troe_first:
                itroe_qss[0] = self.qssReactions.index(reac_id)
                troe_first = False
            itroe_qss[1] = self.qssReactions.index(reac_id) + 1
        if reac_id >= isri[0] and reac_id < isri[1]:
            print(
                "reaction ",
                reac_id,
                mechanism.reaction(id=reac_id).equation(),
                " goes in sri",
            )
            nsri_qss += 1
            if sri_first:
                isri_qss[0] = self.qssReactions.index(reac_id)
                sri_first = False
            isri_qss[1] = self.qssReactions.index(reac_id) + 1
        if reac_id >= ilindemann[0] and reac_id < ilindemann[1]:
            print(
                "reaction ",
                reac_id,
                mechanism.reaction(id=reac_id).equation(),
                " goes in lindemann",
            )
            nlindemann_qss += 1
            if lindemann_first:
                ilindemann_qss[0] = self.qssReactions.index(reac_id)
                lindemann_first = False
            ilindemann_qss[1] = self.qssReactions.index(reac_id) + 1
        if reac_id >= i3body[0] and reac_id < i3body[1]:
            print(
                "reaction ",
                reac_id,
                mechanism.reaction(id=reac_id).equation(),
                " goes in 3body",
            )
            n3body_qss += 1
            if threebody_first:
                i3body_qss[0] = self.qssReactions.index(reac_id)
                threebody_first = False
            i3body_qss[1] = self.qssReactions.index(reac_id) + 1
        if reac_id >= isimple[0] and reac_id < isimple[1]:
            nsimple_qss += 1
            if simple_first:
                isimple_qss[0] = self.qssReactions.index(reac_id)
                simple_first = False
            isimple_qss[1] = self.qssReactions.index(reac_id) + 1
        if reac_id >= ispecial[0] and reac_id < ispecial[1]:
            print(
                "reaction ",
                reac_id,
                mechanism.reaction(id=reac_id).equation(),
                " goes in special",
            )
            nspecial_qss += 1
            if special_first:
                ispecial_qss[0] = self.qssReactions.index(reac_id)
                special_first = False
            ispecial_qss[1] = self.qssReactions.index(reac_id) + 1

    if len(reaction_info.index) != 7:
        print("\n\nCheck this!!!\n")
        sys.exit(1)

    # k_f_qss function
    cw.writer(fstream)
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_k_f_qss(const amrex::Real * tc, amrex::Real invT, amrex::Real * k_f)",
    )
    cw.writer(fstream, "{")
    for index, qssa_reac in enumerate(self.qssReactions):
        reaction = mechanism.reaction()[qssa_reac]
        cw.writer(
            fstream,
            self.line("reaction %d: %s" % (reaction.id, reaction.equation())),
        )
        properties = self.getAllReactionProperties(reaction)
        if properties["fwd_beta"] == 0.0:
            if properties["fwd_Ea"] == 0.0:
                cw.writer(
                    fstream,
                    "k_f[%d] = %.15g * %.15g;"
                    % (
                        index,
                        properties["prefactor_units"],
                        properties["fwd_A"],
                    ),
                )
            else:
                cw.writer(
                    fstream,
                    "k_f[%d] = %.15g * %.15g"
                    % (
                        index,
                        properties["prefactor_units"],
                        properties["fwd_A"],
                    ),
                )
                if properties["activation_units"] > 0:
                    cw.writer(
                        fstream,
                        "           * exp(- %.15g * (%.15g) * invT);"
                        % (
                            properties["activation_units"],
                            properties["fwd_Ea"],
                        ),
                    )
                else:
                    sys.exit(
                        "activation_units neg for reac" + reaction.equation()
                    )
        else:
            cw.writer(
                fstream,
                "k_f[%d] = %.15g * %.15g"
                % (
                    index,
                    properties["prefactor_units"],
                    properties["fwd_A"],
                ),
            )
            if properties["fwd_Ea"] == 0.0:
                cw.writer(
                    fstream,
                    "           * exp(%.15g * tc[0]);"
                    % (properties["fwd_beta"]),
                )
            else:
                if properties["activation_units"] > 0:
                    cw.writer(
                        fstream,
                        "           * exp(%.15g * tc[0] - %.15g * (%.15g) * invT);"
                        % (
                            properties["fwd_beta"],
                            properties["activation_units"],
                            properties["fwd_Ea"],
                        ),
                    )
                else:
                    sys.exit(
                        "activation_units neg for reac" + reaction.equation()
                    )

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")

    # qss coefficients
    cw.writer(fstream)
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qssa_coeff(amrex::Real * k_f, amrex::Real * qf, amrex::Real * qr, amrex::Real * sc, const amrex::Real * tc, amrex::Real * g_RT, amrex::Real * g_RT_qss)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "const amrex::Real invT = 1.0 / tc[1];")
    cw.writer(fstream)

    if len(mechanism.reaction()) == 0:
        cw.writer(fstream)
    else:
        cw.writer(
            fstream,
            self.line(
                "reference concentration: P_atm / (RT) in inverse mol/m^3"
            ),
        )
        cw.writer(
            fstream,
            "const amrex::Real refC = %g / %g * invT;" % (atm.value, R.value),
        )
        cw.writer(fstream, "const amrex::Real refCinv = 1. / refC;")

    cw.writer(fstream, self.line("compute the mixture concentration"))
    cw.writer(fstream, "amrex::Real mixture = 0.0;")
    cw.writer(
        fstream, "for (int i = 0; i < %d; ++i) {" % species_info.n_species
    )
    cw.writer(fstream, "mixture += sc[i];")
    cw.writer(fstream, "}")

    nclassd_qss = self.nqssReactions - nspecial_qss
    nCorr_qss = n3body_qss + ntroe_qss + nsri_qss + nlindemann_qss

    for i in range(nclassd_qss):
        cw.writer(fstream)
        reaction = mechanism.reaction(id=self.qssReactions[i])
        remove_forward = cu.is_remove_forward(reaction_info, orig_idx)
        idx = i
        thirdBody = reaction.thirdBody
        low = reaction.low

        cw.writer(fstream, "{")
        cw.writer(
            fstream,
            self.line("reaction %d: %s" % (reaction.id, reaction.equation())),
        )
        if len(reaction.ford) > 0:
            forward_sc = QSSreturnCoeff(mechanism, reaction.ford)
        else:
            forward_sc = QSSreturnCoeff(mechanism, reaction.reactants)
        if reaction.reversible:
            reverse_sc = QSSreturnCoeff(mechanism, reaction.products)
        else:
            reverse_sc = "0.0"

        KcExpArg = sortedKcExpArg(mechanism, reaction)
        KcConvInv = KcConvInv(mechanism, reaction)

        properties = self.getAllReactionProperties(reaction)

        alpha = 1.0
        if not thirdBody:
            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(
                    fstream,
                    cw.comment(
                        "qf[%d] = k_f[%d] * (%s);" % (idx, idx, forward_sc)
                    ),
                )
                cw.writer(fstream, "qf[%d] = 0.0;" % (idx))
            else:
                cw.writer(
                    fstream,
                    "qf[%d] = k_f[%d] * (%s);" % (idx, idx, forward_sc),
                )
        elif not low:
            alpha = enhancement_d(mechanism, reaction)
            cw.writer(fstream, "const amrex::Real Corr = %s;" % (alpha))
            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(
                    fstream,
                    cw.comment(
                        "qf[%d] = Corr * k_f[%d] * (%s);"
                        % (idx, idx, forward_sc)
                    ),
                )
                cw.writer(fstream, "qf[%d] = 0.0;" % (idx))
            else:
                cw.writer(
                    fstream,
                    "qf[%d] = Corr * k_f[%d] * (%s);" % (idx, idx, forward_sc),
                )
        else:
            alpha = enhancement_d(mechanism, reaction)
            cw.writer(fstream, "amrex::Real Corr = %s;" % (alpha))
            cw.writer(
                fstream,
                "const amrex::Real redP = Corr / k_f[%d] * %.15g "
                % (idx, properties["phase_units"] * properties["low_A"]),
            )
            cw.writer(
                fstream,
                "           * exp(%.15g * tc[0] - %.15g * invT);"
                % (
                    properties["low_beta"],
                    (properties["activation_units"] * properties["low_Ea"]),
                ),
            )
            if reaction.troe:
                cw.writer(
                    fstream, "const amrex::Real F = redP / (1.0 + redP);"
                )
                cw.writer(fstream, "const amrex::Real logPred = log10(redP);")
                cw.writer(fstream, "const amrex::Real logFcent = log10(")
                if abs(properties["troe_Tsss"]) > 1.0e-100:
                    if 1.0 - properties["troe_a"] != 0:
                        cw.writer(
                            fstream,
                            "    %.15g * exp(-tc[1] * %.15g)"
                            % (
                                1.0 - properties["troe_a"],
                                (1 / properties["troe_Tsss"]),
                            ),
                        )
                else:
                    cw.writer(fstream, "     0.0 ")
                if abs(properties["troe_Ts"]) > 1.0e-100:
                    if properties["troe_a"] != 0:
                        cw.writer(
                            fstream,
                            "    + %.15g * exp(-tc[1] * %.15g)"
                            % (
                                properties["troe_a"],
                                (1 / properties["troe_Ts"]),
                            ),
                        )
                else:
                    cw.writer(fstream, "     0.0 ")
                if properties["troe_len"] == 4:
                    if properties["troe_Tss"] < 0:
                        cw.writer(
                            fstream,
                            "    + exp(%.15g * invT));"
                            % -properties["troe_Tss"],
                        )
                    else:
                        cw.writer(
                            fstream,
                            "    + exp(-%.15g * invT));"
                            % properties["troe_Tss"],
                        )
                else:
                    cw.writer(fstream, "    + 0.0);")
                cw.writer(
                    fstream,
                    "const amrex::Real troe_c = -0.4 - 0.67 * logFcent;",
                )
                cw.writer(
                    fstream,
                    "const amrex::Real troe_n = 0.75 - 1.27 * logFcent;",
                )
                cw.writer(
                    fstream,
                    "const amrex::Real troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));",
                )
                cw.writer(
                    fstream,
                    "const amrex::Real F_troe = pow(10, logFcent / (1.0 + troe * troe));",
                )
                cw.writer(fstream, "Corr = F * F_troe;")
                if remove_forward:
                    cw.writer(fstream, cw.comment("Remove forward reaction"))
                    cw.writer(
                        fstream,
                        cw.comment(
                            "qf[%d]  = Corr * k_f[%d] * (%s);"
                            % (idx, idx, forward_sc)
                        ),
                    )
                    cw.writer(fstream, "qf[%d]  = 0.0;" % (idx))
                else:
                    cw.writer(
                        fstream,
                        "qf[%d]  = Corr * k_f[%d] * (%s);"
                        % (idx, idx, forward_sc),
                    )
            elif nlindemann > 0:
                cw.writer(fstream, "Corr = redP / (1.0 + redP);")
                if remove_forward:
                    cw.writer(fstream, cw.comment("Remove forward reaction"))
                    cw.writer(
                        fstream,
                        cw.comment(
                            "qf[%d] = Corr * k_f[%d] * (%s);"
                            % (idx, idx, forward_sc)
                        ),
                    )
                    cw.writer(fstream, "qf[%d] = 0.0;" % (idx))
                else:
                    cw.writer(
                        fstream,
                        "qf[%d] = Corr * k_f[%d] * (%s);"
                        % (idx, idx, forward_sc),
                    )

        if reaction.rev:
            if properties["rev_beta"] == 0:
                if properties["rev_Ea"] == 0:
                    cw.writer(fstream, ";")
                else:
                    cw.writer(
                        fstream,
                        "           * exp( - (%.15g) * invT);"
                        % (
                            properties["activation_units_rev"]
                            * properties["rev_Ea"]
                        ),
                    )
            else:
                if properties["rev_Ea"] == 0:
                    cw.writer(
                        fstream,
                        "           * exp(%.15g * tc[0]);"
                        % (properties["rev_beta"]),
                    )
                else:
                    cw.writer(
                        fstream,
                        "           * exp(%.15g * tc[0] - (%.15g) * invT);"
                        % (
                            properties["rev_beta"],
                            properties["activation_units_rev"]
                            * properties["rev_Ea"],
                        ),
                    )

            if alpha == 1.0:
                cw.writer(
                    fstream,
                    "qr[%d] = k_r[%d] * (%s);" % (idx, idx, reverse_sc),
                )
            else:
                cw.writer(
                    fstream,
                    "qr[%d] = Corr * k_r[%d] * (%s);" % (idx, idx, reverse_sc),
                )
        else:
            if KcConvInv:
                if alpha == 1.0:
                    cw.writer(
                        fstream,
                        "qr[%d] = k_f[%d] * exp(-(%s)) * (%s) * (%s);"
                        % (idx, idx, KcExpArg, KcConvInv, reverse_sc),
                    )
                else:
                    cw.writer(
                        fstream,
                        "qr[%d] = Corr * k_f[%d] * exp(-(%s)) * (%s) * (%s);"
                        % (idx, idx, KcExpArg, KcConvInv, reverse_sc),
                    )
            else:
                if alpha == 1.0:
                    cw.writer(
                        fstream,
                        "qr[%d] = k_f[%d] * exp(-(%s)) * (%s);"
                        % (idx, idx, KcExpArg, reverse_sc),
                    )
                else:
                    cw.writer(
                        fstream,
                        "qr[%d] = Corr * k_f[%d] * exp(-(%s)) * (%s);"
                        % (idx, idx, KcExpArg, reverse_sc),
                    )
        cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")

    # qss concentrations
    cw.writer(fstream)
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_sc_qss(amrex::Real * sc_qss, amrex::Real * qf_co, amrex::Real * qr_co)",
    )
    cw.writer(fstream, "{")

    cw.writer(fstream, "amrex::Real epsilon = 1e-12;")
    cw.writer(fstream)

    # print()
    print("** self.decouple_index:")
    print(self.decouple_index)
    print(list(self.needs.keys()))
    print(list(self.group.keys()))

    for i in self.decouple_index:
        symbol = self.decouple_index[i]
        print("... Dealing with Spec/Group ", symbol)
        if symbol in list(self.needs.keys()):
            print("    Simple case, single group")

            denominator = symbol + "_denom"
            numerator = symbol + "_num"

            cw.writer(
                fstream,
                self.line(
                    "QSSA species "
                    + str(species_info.qssa_species_list.index(symbol))
                    + ": "
                    + symbol
                ),
            )
            cw.writer(fstream)
            # RHS
            # cut line if too big !
            long_line_elements = (self.qssa_rhs[symbol]).split()
            len_long_line = len(long_line_elements)
            # if we have more than 7 elements
            if len_long_line > 7:
                # treat first line separately with the epsilon
                cw.writer(
                    fstream,
                    "amrex::Real %s = epsilon %s"
                    % (numerator, " ".join(long_line_elements[0:7])),
                )
                # proceed by strides of 7
                for kk in range(7, len_long_line, 7):
                    # if there are less than 7 elems left then we are at the end of the list
                    if len(long_line_elements[kk : kk + 7]) < 7:
                        cw.writer(
                            fstream,
                            "                    %s;"
                            % (" ".join(long_line_elements[kk : kk + 7])),
                        )
                    # if there are 7 elems we are ...
                    else:
                        # either are in the middle of the list
                        if len(long_line_elements[kk:]) > 7:
                            cw.writer(
                                fstream,
                                "                    %s"
                                % (" ".join(long_line_elements[kk : kk + 7])),
                            )
                        # or at the end but list number was a multiple of 7
                        else:
                            cw.writer(
                                fstream,
                                "                    %s;"
                                % (" ".join(long_line_elements[kk : kk + 7])),
                            )
            # if we have less than 7 elements just write them
            else:
                cw.writer(
                    fstream,
                    "amrex::Real %s = epsilon %s;"
                    % (numerator, self.qssa_rhs[symbol]),
                )
            # COEFF
            cw.writer(
                fstream,
                "amrex::Real %s = epsilon %s;"
                % (denominator, self.qssa_coeff[symbol]),
            )
            cw.writer(fstream)
            cw.writer(
                fstream,
                "sc_qss[%s] = - %s/%s;"
                % (
                    species_info.qssa_species_list.index(symbol),
                    numerator,
                    denominator,
                ),
            )
            cw.writer(fstream)
        if symbol in list(self.group.keys()):
            print(
                "    Though case. Submatrix has size ",
                len(self.group[symbol]),
                "x",
                len(self.group[symbol]),
            )
            coeff_submatrix = [
                ["0"] * len(self.group[symbol])
                for i in range(len(self.group[symbol]))
            ]
            rhs_submatrix = ["0"] * len(self.group[symbol])
            gr_species = self.group[symbol]
            print("    Species involved :", gr_species)
            cw.writer(
                fstream,
                "/* QSSA coupling between " + ("  ").join(gr_species) + "*/",
            )
            for index, species in enumerate(gr_species):
                print("      x Dealing with spec", species, index)
                cw.writer(
                    fstream,
                    self.line(
                        "QSSA species "
                        + str(species_info.qssa_species_list.index(species))
                        + ": "
                        + species
                    ),
                )
                cw.writer(fstream)

                denominator = species + "_denom"
                numerator = species + "_num"

                # RHS
                # cut line if too big !
                long_line_elements = (self.qssa_rhs[species]).split()
                len_long_line = len(long_line_elements)
                # if we have more than 7 elements
                if len_long_line > 7:
                    # treat first line separately with the epsilon
                    cw.writer(
                        fstream,
                        "amrex::Real %s = epsilon %s"
                        % (numerator, " ".join(long_line_elements[0:7])),
                    )
                    # proceed by strides of 7
                    for kk in range(7, len_long_line, 7):
                        # if there are less than 7 elems left then we are at the end of the list
                        if len(long_line_elements[kk : kk + 7]) < 7:
                            cw.writer(
                                fstream,
                                "                    %s;"
                                % (" ".join(long_line_elements[kk : kk + 7])),
                            )
                        # if there are 7 elems we are ...
                        else:
                            # either are in the middle of the list
                            if len(long_line_elements[kk:]) > 7:
                                cw.writer(
                                    fstream,
                                    "                    %s"
                                    % (
                                        " ".join(
                                            long_line_elements[kk : kk + 7]
                                        )
                                    ),
                                )
                            # or at the end but list number was a multiple of 7
                            else:
                                cw.writer(
                                    fstream,
                                    "                    %s;"
                                    % (
                                        " ".join(
                                            long_line_elements[kk : kk + 7]
                                        )
                                    ),
                                )
                # if we have less than 7 elements just write them
                else:
                    cw.writer(
                        fstream,
                        "amrex::Real %s = epsilon %s;"
                        % (numerator, self.qssa_rhs[species]),
                    )
                # COEFF
                # cut line if too big !
                long_line_elements = (self.qssa_coeff[species]).split()
                len_long_line = len(long_line_elements)
                # if we have more than 7 elements
                if len_long_line > 7:
                    # treat first line separately with the epsilon
                    cw.writer(
                        fstream,
                        "amrex::Real %s = epsilon %s"
                        % (denominator, " ".join(long_line_elements[0:7])),
                    )
                    # proceed by strides of 7
                    for kk in range(7, len_long_line, 7):
                        # if there are less than 7 elems left then we are at the end of the list
                        if len(long_line_elements[kk : kk + 7]) < 7:
                            cw.writer(
                                fstream,
                                "                    %s;"
                                % (" ".join(long_line_elements[kk : kk + 7])),
                            )
                        # if there are 7 elems we are ...
                        else:
                            # either are in the middle of the list
                            if len(long_line_elements[kk:]) > 7:
                                cw.writer(
                                    fstream,
                                    "                    %s"
                                    % (
                                        " ".join(
                                            long_line_elements[kk : kk + 7]
                                        )
                                    ),
                                )
                            # or at the end but list number was a multiple of 7
                            else:
                                cw.writer(
                                    fstream,
                                    "                    %s;"
                                    % (
                                        " ".join(
                                            long_line_elements[kk : kk + 7]
                                        )
                                    ),
                                )
                # if we have less than 7 elements just write them
                else:
                    cw.writer(
                        fstream,
                        "amrex::Real %s = epsilon %s;"
                        % (denominator, self.qssa_coeff[species]),
                    )
                # RHS
                cw.writer(
                    fstream,
                    "amrex::Real "
                    + species
                    + "_rhs = -"
                    + numerator
                    + "/"
                    + denominator
                    + ";",
                )
                cw.writer(fstream)

                for j in range(len(gr_species)):
                    if j == index:
                        coeff_submatrix[index][j] = "1"
                    else:
                        if (
                            self.qssa_qssa_coeff[species][gr_species[j]]
                            != "0.0"
                        ):
                            coeff_submatrix[index][j] = (
                                str(species) + "_" + str(gr_species[j])
                            )
                            # let us assume for now these lines are not too big
                            cw.writer(
                                fstream,
                                "amrex::Real "
                                + str(species)
                                + "_"
                                + str(gr_species[j])
                                + " = (epsilon "
                                + self.qssa_qssa_coeff[species][gr_species[j]]
                                + ")/"
                                + denominator
                                + ";",
                            )
                cw.writer(fstream)
                rhs_submatrix[index] = str(species) + "_rhs"

            A, X, B, intermediate_helpers = gauss_pivoting(
                species_info, coeff_submatrix, rhs_submatrix
            )

            print("X is ", X)

            cw.writer(fstream, cw.comment("Putting it all together"))
            for helper in intermediate_helpers:
                if helper in self.list_of_intermediate_helpers:
                    cw.writer(
                        fstream,
                        "%s = %s;" % (helper, intermediate_helpers[helper]),
                    )
                else:
                    cw.writer(
                        fstream,
                        "amrex::Real %s = %s;"
                        % (helper, intermediate_helpers[helper]),
                    )
                    self.list_of_intermediate_helpers.append(helper)

            for count in range(len(gr_species)):
                max_index = len(gr_species) - 1
                species = gr_species[max_index - count]

                # cut line if too big !
                long_line_elements = X[max_index - count].split()
                len_long_line = len(long_line_elements)
                # if we have more than 4 elements
                if len_long_line > 4:
                    # treat first line separately
                    cw.writer(
                        fstream,
                        "sc_qss["
                        + str(species_info.qssa_species_list.index(species))
                        + "] = "
                        + (" ".join(long_line_elements[0:4])),
                    )
                    # proceed by strides of 4
                    for kk in range(4, len_long_line, 4):
                        # if there are less than 4 elems left then we are at the end of the list
                        if len(long_line_elements[kk : kk + 4]) < 4:
                            cw.writer(
                                fstream,
                                "                    %s;"
                                % (" ".join(long_line_elements[kk : kk + 4])),
                            )
                        # if there are 4 elems we are ...
                        else:
                            # either are in the middle of the list
                            if len(long_line_elements[kk:]) > 4:
                                cw.writer(
                                    fstream,
                                    "                    %s"
                                    % (
                                        " ".join(
                                            long_line_elements[kk : kk + 4]
                                        )
                                    ),
                                )
                            # or at the end but list number was a multiple of 4
                            else:
                                cw.writer(
                                    fstream,
                                    "                    %s;"
                                    % (
                                        " ".join(
                                            long_line_elements[kk : kk + 4]
                                        )
                                    ),
                                )
                # if we have less than 4 elements just write them
                else:
                    cw.writer(
                        fstream,
                        "sc_qss["
                        + str(species_info.qssa_species_list.index(species))
                        + "] = "
                        + X[max_index - count]
                        + ";",
                    )
                cw.writer(fstream)

    cw.writer(fstream)
    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")

    return


def gauss_pivoting(species_info, A, B=None):
    """Gauss pivoting."""
    print()
    print("In Gauss pivot")

    pivots = []
    intermediate_helpers = OrderedDict()
    helper_counters = 0

    X = [""] * len(A[0])
    for i in range(len(A[0])):
        X[i] = "X" + str(i)

    if B is None:
        for i in range(len(A[0])):
            B[i] = "B" + str(i)

    # Get species names:
    species = ["0" for i in range(len(B))]
    for member in range(len(B)):
        hold = str(B[member])
        hold = hold[:-4]
        species[member] = hold

    print("Species involved are: ", species)
    print()

    Anum = np.zeros([len(A[0]), len(A[0])])
    for i in range(len(A[0])):
        for j in range(len(A[0])):
            if A[i][j] != "0":
                Anum[i, j] = 1

    print("--A", A)
    print("--B", B)

    indi, indj = np.nonzero(Anum)

    n = len(B)
    for k in range(n - 1):

        pivot = A[k][k]

        # swap lines if needed
        if pivot == 0:
            temp = np.array(A[k + 1][:])
            A[k + 1][:] = A[k][:]
            A[k][:] = temp

            temp = str(B[k + 1])
            B[k + 1] = B[k]
            B[k] = temp

            pivot = A[k][k]

        print()
        print("   **ROW of pivot ", k, " and pivot is ", pivot)
        pivots.append(pivot)

        for i in range(k, len(B) - 1):
            num = A[i + 1][k]
            print(
                "       xx Treating ROW ",
                i + 1,
                "with numerator ",
                num,
                "subst fact will be",
                num,
                "/",
                pivot,
            )
            if num == "0":
                print(
                    "        !! No need to do anything, already zeroed ... skip"
                )
                continue
            B = list(B)
            print("          - B starts with: ")
            print("           ", B)
            if num != "0":
                if pivot != "1":
                    if num != "1":
                        helper = num + "/" + pivot
                        helper_name = "H_" + str(helper_counters)
                        intermediate_helpers[helper_name] = helper
                        B[i + 1] = (
                            B[i + 1] + " -" + B[int(k)] + "*" + helper_name
                        )
                        B[i + 1] = "(" + B[i + 1] + ")"
                        helper_counters += 1
                    else:
                        helper = 1 + "/" + pivot
                        helper_name = "H_" + str(helper_counters)
                        intermediate_helpers[helper_name] = helper
                        print(" IN THIS CASE !! CHECK THAT ITS OK !! ")
                        B[i + 1] = (
                            B[i + 1] + " -" + B[int(k)] + "*" + helper_name
                        )
                        B[i + 1] = "(" + B[i + 1] + ")"
                        helper_counters += 1
                else:
                    if num != "1":
                        helper = num
                        helper_name = "H_" + str(helper_counters)
                        intermediate_helpers[helper_name] = helper
                        B[i + 1] = (
                            B[i + 1] + " -" + B[int(k)] + "*" + helper_name
                        )
                        B[i + 1] = "(" + B[i + 1] + ")"
                        helper_counters += 1
                    else:
                        B[i + 1] = B[i + 1] + " -" + B[int(k)]
                        B[i + 1] = "(" + B[i + 1] + ")"

            print("          ... and B ends with: ")
            print("            ", B)

            indi, indj = np.nonzero(Anum)

            for j in indj[indi == (i + 1)]:
                print(
                    "          - Dealing with row elem on column ",
                    j,
                    " : ",
                    A[i + 1][j],
                )
                if j == k:
                    print("            !! 0 ELEM !")
                    A[i + 1][j] = "0"
                else:
                    if A[i + 1][j] != "0":
                        if num != "0":
                            if pivot != "1":
                                if num != "1":
                                    if A[k][j] != "0":
                                        A[i + 1][j] = (
                                            A[i + 1][j]
                                            + " -"
                                            + A[k][j]
                                            + "*"
                                            + helper_name
                                        )
                                        A[i + 1][j] = "(" + A[i + 1][j] + ")"
                                else:
                                    if A[k][j] != "0":
                                        print(
                                            " IN THIS CASE !! CHECK THAT ITS OK !! "
                                        )
                                        A[i + 1][j] = (
                                            A[i + 1][j]
                                            + " -"
                                            + A[k][j]
                                            + "*"
                                            + helper_name
                                        )
                                        A[i + 1][j] = "(" + A[i + 1][j] + ")"
                            else:
                                if num != "1":
                                    if A[k][j] != "0":
                                        A[i + 1][j] = (
                                            A[i + 1][j]
                                            + " -"
                                            + A[k][j]
                                            + "*"
                                            + helper_name
                                        )
                                        A[i + 1][j] = "(" + A[i + 1][j] + ")"
                                else:
                                    if A[k][j] != "0":
                                        A[i + 1][j] = (
                                            A[i + 1][j] + " -" + A[k][j]
                                        )
                                        A[i + 1][j] = "(" + A[i + 1][j] + ")"
                    else:
                        if num != "0":
                            if pivot != "1":
                                if num != "1":
                                    if A[k][j] != "0":
                                        A[i + 1][j] = (
                                            " -" + A[k][j] + "*" + helper_name
                                        )
                                        A[i + 1][j] = "(" + A[i + 1][j] + ")"
                                else:
                                    if A[k][j] != "0":
                                        print(
                                            " IN THIS CASE !! CHECK THAT ITS OK !! "
                                        )
                                        A[i + 1][j] = (
                                            " -" + A[k][j] + "*" + helper_name
                                        )
                                        A[i + 1][j] = "(" + A[i + 1][j] + ")"
                            else:
                                if num != "1":
                                    if A[k][j] != "0":
                                        A[i + 1][j] = (
                                            " -" + A[k][j] + "*" + helper_name
                                        )
                                        A[i + 1][j] = "(" + A[i + 1][j] + ")"
                                else:
                                    if A[k][j] != "0":
                                        A[i + 1][j] = " -" + A[k][j]
                                        A[i + 1][j] = "(" + A[i + 1][j] + ")"
            print("          ... and updated A is: ")
            print("             ", A)

    for i in range(len(B)):
        X = list(X)
        B[i] = str(B[i])

    # start with last elem
    n = n - 1
    if A[n][n] != "1":
        X[n] = B[n] + "/" + A[n][n]
    else:
        X[n] = B[n]

    for i in range(1, n + 1):
        sumprod = ""
        for j in range(i):
            flag = False
            if A[n - i][n - j] != "0":
                if flag:
                    sumprod += " + "
                flag = True
                if A[n - i][n - j] == "1":
                    sumprod += " (" + str(X[n - j]) + ")"
                elif j != 0:
                    sumprod += (
                        " +"
                        + A[n - i][n - j]
                        + "*"
                        + "sc_qss["
                        + str(
                            species_info.qssa_species_list.index(
                                species[n - j]
                            )
                        )
                        + "]"
                    )
                else:
                    sumprod += (
                        A[n - i][n - j]
                        + "*"
                        + "sc_qss["
                        + str(
                            species_info.qssa_species_list.index(
                                species[n - j]
                            )
                        )
                        + "]"
                    )

        if sumprod == "":
            if A[n - i][n - i] != "1":
                X[n - i] = "(" + B[n - i] + ")/" + A[n - i][n - i]
            else:
                X[n - i] = B[n - i]
        else:
            if A[n - i][n - i] == "1":
                X[n - i] = B[n - i] + " -(" + sumprod + ")"
            else:
                X[n - i] = (
                    "(" + B[n - i] + " -(" + sumprod + "))/" + A[n - i][n - i]
                )
    print()
    print()

    return A, X, B, intermediate_helpers
