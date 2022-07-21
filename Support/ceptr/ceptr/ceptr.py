"""Convert cantera mechanism to C++ files."""
import argparse
from argparse import RawTextHelpFormatter
import pathlib

import cantera as ct

import ceptr.converter as converter


def convert(
    fname,
    hformat,
    remove_1,
    remove_pow,
    remove_pow10,
    min_op_count,
    gradual_op_count,
    store_in_jacobian,
    round_decimals,
    recycle_cse,
    min_op_count_all,
    remove_single_symbols_cse,
):
    """Convert a mechanism file."""
    mechanism = ct.Solution(fname)
    conv = converter.Converter(
        mechanism,
        hformat,
        remove_1,
        remove_pow,
        remove_pow10,
        min_op_count,
        gradual_op_count,
        store_in_jacobian,
        round_decimals,
        recycle_cse,
        min_op_count_all,
        remove_single_symbols_cse,
    )
    conv.writer()
    conv.formatter()


def convert_lst(lst, hformat):
    """Convert mechanisms from a file containing a list of directories."""
    lpath = pathlib.Path(lst)
    with open(lst, "r") as f:
        for line in f:
            if not line.startswith("#"):
                mechname = lpath.parents[0] / line.strip()
                print(f"""Converting file {mechname}""")
                convert(mechname, hformat)


def main():
    """Convert cantera mechanisms to C++ files."""
    parser = argparse.ArgumentParser(description="Mechanism converter", formatter_class=RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fname", help="Mechanism file", type=str)
    group.add_argument(
        "-l", "--lst", help="Mechanism directory file list", type=str
    )

    long_help = "Sytle format for .H file output."
    long_help += "\nCPU: will print intermediate variables used for chainruling. This gives a readable version of the Jacobian entries, albeit memory consuming."
    long_help += "\nGPU: will not print intermediate variables used for chainruling, and instead will replace them directly in the Jacobian entries."
    long_help += " This gives a less readable version of the Jacobian, but more memory efficient."
    parser.add_argument(
        "--hformat",
        help=long_help,
        type=str,
        choices=["cpu", "gpu"],
        default="gpu",
        required=False,
    )

    parser.add_argument(
        "-r1",
        "--remove_1",
        action="store_true",
        help="Remove factor 1.0 in printed expressions",
    )

    parser.add_argument(
        "-rp",
        "--remove_pow",
        action="store_true",
        help="Replace pow(...,n) with multiplications or divisions if n<=3 and n>=-3 in printed expressions.",
    )

    parser.add_argument(
        "-rp10",
        "--remove_pow10",
        action="store_true",
        help="Remove pow(10,x) in printed expressions and replace it with exp(ln(10)*x).",
    )

    long_help = "Counts number operations used to construct each common subexpression and replace"
    long_help += " the common subexpression if the number of operations is less or equal to the value"
    parser.add_argument(
        "-moc",
        "--min_op_count",
        type=int,
        metavar="",
        required=False,
        help=long_help,
        default=0,
    )

    long_help = "Similar to --min_op_count but also counts how many times that common subexpression is used later."
    long_help += "\nThe meaning value passed is how many more operations will be done if the common subexpression is eliminated."
    long_help += "\nThis option only marginally increase the file size (therefore compile time), while still being memory efficient."
    parser.add_argument(
        "-moca",
        "--min_op_count_all",
        type=int,
        metavar="",
        required=False,
        help=long_help,
        default=0,
    )

    long_help = "Gradual elimination of common subexpressions."
    long_help += "\nUseful if --min_op_count or --min_op_count_all are active."
    long_help += "\nLoops from 1 to the min_op_count and min_op_count_all values and gradually eliminate the common subexpressions."
    long_help += "\nThis has the advantage of ensuring that the memory footprint is strictly monotonically decreasing as min_op_count and min_op_count_all are increased."
    parser.add_argument(
        "-roc",
        "--gradual_op_count",
        action="store_true",
        help=long_help,
    )

    long_help = "Use the Jacobian array as a temporary space to store intermediate variables."
    long_help += "\nIn particular, the last row of the Jacobian (dependence with respect to temperature) is done by finite difference which requires storing intermediate variables"
    long_help += " (production rate, forward and backward reactions)."
    long_help += "\nWhen the option is active, the `productionRate` function used to compute the finite difference"
    long_help += " is replaced with a `productionRate_light` functions where references to different parts of the Jacobian are used in place of allocating new arrays."
    parser.add_argument(
        "-sj",
        "--store_in_jacobian",
        action="store_true",
        help=long_help,
    )

    parser.add_argument(
        "-rd",
        "--round_decimals",
        action="store_true",
        help="Round decimal numbers when possible to minimize character count",
    )

    parser.add_argument(
        "-rcse",
        "--recycle_cse",
        action="store_true",
        help="Reuse common subexpressions that are not used later to avoid declaring new temporary reals",
    )

    long_help = "Remove common subexpressions that are made of 1 symbol."
    long_help += "\nThose common subexpressions are typically `-xxx` and may not appear as worth replacing because they save 1 operations and are reused multiple times."
    long_help += "\nHowever, when replaced in the later expressions, the `-` operations typically disappear or is merged into another"
    long_help += " operations which actually does not increase the total number of operations."
    parser.add_argument(
        "-rss",
        "--remove_single_symbols_cse",
        action="store_true",
        help=long_help,
    )

    args = parser.parse_args()

    if args.fname:
        convert(
            args.fname,
            args.hformat,
            args.remove_1,
            args.remove_pow,
            args.remove_pow10,
            args.min_op_count,
            args.gradual_op_count,
            args.store_in_jacobian,
            args.round_decimals,
            args.recycle_cse,
            args.min_op_count_all,
            args.remove_single_symbols_cse,
        )
    elif args.lst:
        convert_lst(
            args.lst,
            args.hformat,
            args.remove_1,
            args.remove_pow,
            args.remove_pow10,
            args.min_op_count,
            args.gradual_op_count,
            args.store_in_jacobian,
            args.round_decimals,
            args.recycle_cse,
            args.min_op_count_all,
            args.remove_single_symbols_cse,
        )


if __name__ == "__main__":
    main()
