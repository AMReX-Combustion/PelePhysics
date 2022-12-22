"""Convert cantera mechanism to C++ files."""
import argparse
import pathlib

import cantera as ct

import ceptr.converter as converter


def convert(
    fname,
    jocobian,
    qss_format_input,
    qss_symbolic_jac,
):
    """Convert a mechanism file."""
    mechanism = ct.Solution(fname)
    conv = converter.Converter(
        mechanism,
        jocobian,
        qss_format_input,
        qss_symbolic_jac,
    )
    conv.writer()
    conv.formatter()


def convert_lst(
    lst,
    jocobian,
    qss_format_input,
    qss_symbolic_jac,
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
                    jocobian,
                    qss_format_input,
                    qss_symbolic_jac,
                )


def main():
    """Convert cantera mechanisms to C++ files."""
    parser = argparse.ArgumentParser(
        description="Mechanism converter",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fname", help="Mechanism file", type=str)
    group.add_argument(
        "-l", "--lst", help="Mechanism directory file list", type=str
    )

    parser.add_argument(
        "--qss_format_input",
        help="Input file for QSS Jacobian formatting parameters mechanisms",
        type=str,
        default=None,
        required=False,
    )

    parser.add_argument(
        "-qsj",
        "--qss_symbolic_jacobian",
        action="store_true",
        help="Compute the QSS Jacobian using symbolic recording",
    )

    parser.add_argument(
        "-nj",
        "--no_jocobian",
        action="store_true",
        help="Do not generate a jacobian",
    )

    args = parser.parse_args()

    if args.fname:
        convert(
            args.fname,
            not args.no_jocobian,
            args.qss_format_input,
            args.qss_symbolic_jacobian,
        )
    elif args.lst:
        convert_lst(
            args.lst,
            not args.no_jocobian,
            args.qss_format_input,
            args.qss_symbolic_jacobian,
        )


if __name__ == "__main__":
    main()
