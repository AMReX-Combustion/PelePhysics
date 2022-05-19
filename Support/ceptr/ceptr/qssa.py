"""Generate QSSA chemistry file."""
import argparse
import pathlib
import sys

import cantera as ct
import yaml

import ceptr.qssa_reduction as cqr
import ceptr.reaction_info as cri


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
    parser.add_argument(
        "-v",
        "--visualize",
        help="Visualize quadratic coupling and QSSA dependencies",
        action="store_true",
    )
    args = parser.parse_args()

    # Load mechanism
    mechanism = ct.Solution(args.fname)
    reaction_info = cri.sort_reactions(mechanism)
    mechpath = pathlib.Path(mechanism.source)
    qssaname = mechpath.parents[0] / "qssa.yaml"

    # Species
    with open(args.nqssa) as f:
        f_non_qssa_species = yaml.safe_load(f)
        non_qssa_species = f_non_qssa_species["species"]
    all_species = mechanism.species_names
    qssa_species = list(set(all_species) - set(non_qssa_species))

    # Visualize
    if args.visualize:
        cqr.visualize_qssa(mechanism, reaction_info, qssa_species)

    forward_to_remove = []
    if args.method == 0:
        qssa_species = cqr.remove_quadratic_method_0(mechanism, qssa_species)
        reactions_to_keep = mechanism.reactions()
    elif args.method == 1 or args.method == 2:

        candidates_for_removal = cqr.identify_removals(
            mechanism, reaction_info, qssa_species
        )

        if args.method == 1:
            reactions_to_keep = cqr.remove_quadratic_method_1(
                mechanism, reaction_info, candidates_for_removal
            )
        elif args.method == 2:
            (
                reactions_to_keep,
                forward_to_remove,
            ) = cqr.remove_quadratic_method_2(
                mechanism, reaction_info, candidates_for_removal
            )

    qssa = ct.Solution(
        name=f"{mechanism.name}-qssa",
        thermo=mechanism.thermo_model,
        kinetics=mechanism.kinetics_model,
        transport_model=mechanism.transport_model,
        species=mechanism.species(),
        reactions=reactions_to_keep,
    )

    if cqr.qssa_coupling(qssa, qssa_species, forward_to_remove):
        print("Failure to generate QSSA mechanism.")
        sys.exit(1)

    qssa.update_user_header({"description": f"QSSA of {mechanism.name}"})
    qssa.update_user_data(
        {"qssa_species": qssa_species, "n_qssa_species": len(qssa_species)}
    )
    forward_to_remove_idx = []
    if forward_to_remove:
        for fr in forward_to_remove:
            for idx, reaction in enumerate(qssa.reactions()):
                if fr.equation == reaction.equation:
                    forward_to_remove_idx.append(idx)
                    break
        qssa.update_user_data({"forward_to_remove_idx": forward_to_remove_idx})
    qssa.write_yaml(qssaname, header=True)

    # Summarize
    for idx in forward_to_remove_idx:
        print(f"Forward reaction to be removed: {qssa.reaction(idx)}")

    mechanism = ct.Solution(args.fname)  # reread
    for reaction in mechanism.reactions():
        if reaction.equation not in [x.equation for x in qssa.reactions()]:
            # check if the unreversed reaction
            reaction.reversible = False
            if reaction.equation in [x.equation for x in qssa.reactions()]:
                reaction.reversible = True
                print("Removed reverse reaction:", reaction)
            else:
                reaction.reversible = True
                print("Removed reaction:", reaction)


if __name__ == "__main__":
    main()
