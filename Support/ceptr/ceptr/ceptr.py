"""Convert cantera mechanism to C++ files."""
import argparse
import pathlib
from argparse import RawTextHelpFormatter

import cantera as ct

import ceptr.converter as converter


def convert(
    fname,
    format_input,
    symbolic_jac,
):
    """Convert a mechanism file."""
    mechanism = ct.Solution(fname)
    conv = converter.Converter(
        mechanism,
        format_input,
        symbolic_jac,
    )
    conv.writer()
    conv.formatter()


def convert_lst(
    lst,
    format_input,
    symbolic_jac,
):
    """Convert mechanisms from a file containing a list of directories."""
    lpath = pathlib.Path(lst)
    with open(lst, "r") as f:
        for line in f:
            if not line.startswith("#"):
                mechname = lpath.parents[0] / line.strip()
                print(f"""Converting file {mechname}""")
                convert(
                    mechname,
                    format_input,
                    symbolic_jac,
                )


def main():
    """Convert cantera mechanisms to C++ files."""
    parser = argparse.ArgumentParser(
        description="Mechanism converter", formatter_class=RawTextHelpFormatter
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fname", help="Mechanism file", type=str)
    group.add_argument(
        "-l", "--lst", help="Mechanism directory file list", type=str
    )

    parser.add_argument(
        "--format_input",
        help=(
            "Input file for jacobian formatting parameters used with QSS"
            " mechanisms"
        ),
        type=str,
        default=None,
        required=False,
    )

    parser.add_argument(
        "-sj",
        "--symbolic_jacobian",
        action="store_true",
        help="Compute the Jacobian using symbolic recording",
    )

    args = parser.parse_args()

    if args.fname:
        convert(
            args.fname,
            args.format_input,
            args.symbolic_jacobian,
        )
    elif args.lst:
        convert_lst(
            args.lst,
            args.format_input,
            args.symbolic_jacobian,
        )


if __name__ == "__main__":
    main()
