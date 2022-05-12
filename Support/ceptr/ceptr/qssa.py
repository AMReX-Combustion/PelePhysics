"""Convert cantera mechanism to C++ files."""
import argparse
import pathlib

import cantera as ct
import yaml

import ceptr.reaction_info as cri


def intersection(lst1, lst2):
    """Return intersection of two lists."""
    return list(set(lst1).intersection(lst2))


def identify_removals(mechanism, reaction_info, qssa_species):
    """Identify the reactions that should be removed."""
    reactions_to_remove = {}
    for orig_idx, _ in reaction_info.idxmap.items():
        reaction = mechanism.reaction(orig_idx)

        # Check for quadratic relations
        qssa_species_involved = intersection(
            list(reaction.reactants.keys()), qssa_species
        )
        sum_coeff = sum(
            v
            for k, v in reaction.reactants.items()
            if k in qssa_species_involved
        )
        if sum_coeff > 1:
            if orig_idx not in reactions_to_remove:
                reactions_to_remove[orig_idx] = []
            reactions_to_remove[orig_idx].append("f")
        if reaction.reversible:
            qssa_species_involved = intersection(
                list(reaction.products.keys()), qssa_species
            )
            sum_coeff = sum(
                v
                for k, v in reaction.products.items()
                if k in qssa_species_involved
            )
            if sum_coeff > 1:
                if orig_idx not in reactions_to_remove:
                    reactions_to_remove[orig_idx] = []
                reactions_to_remove[orig_idx].append("r")

    return reactions_to_remove


def main():
    """Apply QSSA reduction to Cantera yaml file."""
    parser = argparse.ArgumentParser(description="Mechanism converter")
    parser.add_argument(
        "-f", "--fname", help="Mechanism file", type=str, required=True
    )
    parser.add_argument(
        "-n", "--nqssa", help="Non-QSSA species list", type=str, required=True
    )
    parser.add_argument(
        "-m",
        "--method",
        help="QSSA method (default: 2)",
        type=int,
        choices=[0, 1, 2],
        default=2,
    )
    args = parser.parse_args()

    # Load mechanism
    mechanism = ct.Solution(args.fname)
    reaction_info = cri.sort_reactions(mechanism)
    mechpath = pathlib.Path(mechanism.source)
    qssaname = mechpath.parents[0] / "qssa.yaml"

    # Species
    with open(args.nqssa) as f:
        f_species_non_qssa = yaml.safe_load(f)
        species_non_qssa = f_species_non_qssa["species"]
    species_all = mechanism.species_names
    species_qssa = list(set(species_all) - set(species_non_qssa))

    reactions_to_remove = identify_removals(
        mechanism, reaction_info, species_qssa
    )

    print("Removing the following reactions:")
    for k in reactions_to_remove.keys():
        print(mechanism.reaction(k))

    reactions_to_keep = []
    for orig_idx, _ in reaction_info.idxmap.items():
        reaction = mechanism.reaction(orig_idx)
        if orig_idx in reactions_to_remove:
            if (
                "f" in reactions_to_remove[orig_idx]
                and "r" in reactions_to_remove[orig_idx]
            ):
                continue
            elif "f" in reactions_to_remove[orig_idx]:
                if reaction.reversible:
                    if (
                        not hasattr(reaction, "low_rate")
                        and not reaction.reaction_type == "three-body"
                    ):
                        continue
                    else:
                        continue
                else:
                    continue
            elif "r" in reactions_to_remove[orig_idx]:
                reaction.reversible = False
                reactions_to_keep.append(reaction)
        else:
            reactions_to_keep.append(reaction)

    qssa = ct.Solution(
        name=f"{mechanism.name}-qssa",
        thermo=mechanism.thermo_model,
        kinetics=mechanism.kinetics_model,
        transport_model=mechanism.transport_model,
        species=mechanism.species(),
        reactions=reactions_to_keep,
    )
    qssa.update_user_header({"description": f"QSSA of {mechanism.name}"})
    qssa.write_yaml(qssaname, header=True)


if __name__ == "__main__":
    main()
