"""Convert cantera mechanism to C++ files."""
import argparse
import pathlib
import time

import cantera as ct

import ceptr.converter as converter


def convert(
    fname,
    jacobian,
    qss_format_input,
    qss_symbolic_jac,
):
    """Convert a mechanism file."""
    mechanism = ct.Solution(fname)
    conv = converter.Converter(
        mechanism,
        jacobian,
        qss_format_input,
        qss_symbolic_jac,
    )
    conv.writer()
    conv.formatter()


def convert_lst(
    lst,
    jacobian,
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
                    jacobian,
                    qss_format_input,
                    qss_symbolic_jac,
                )


def convert_lst_qss(
    lst,
    jacobian,
):
    """Convert QSS mechanisms from a file of directories and format input."""
    lpath = pathlib.Path(lst)
    with open(lst, "r") as f:
        for line in f:
            if not line.startswith("#"):
                mech_file, format_file, _, _ = line.split()
                mechname = lpath.parents[0] / mech_file.strip()
                qss_format_input = lpath.parents[0] / format_file.strip()
                print(f"""Converting file {mechname}""")
                convert(
                    mechname,
                    jacobian,
                    qss_format_input,
                    True,
                )


def main():
    """Convert cantera mechanisms to C++ files."""
    start = time.time()
    parser = argparse.ArgumentParser(
        description="Mechanism converter",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fname", help="Mechanism file", type=str)
    group.add_argument(
        "-l", "--lst", help="Mechanism directory file list", type=str
    )
    group.add_argument(
        "-lq", "--lst_qss", help="QSS mechanism directory file list", type=str
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
        "--no_jacobian",
        action="store_true",
        help="Do not generate a jacobian",
    )

    args = parser.parse_args()

    if args.fname:
        convert(
            args.fname,
            not args.no_jacobian,
            args.qss_format_input,
            args.qss_symbolic_jacobian,
        )
    elif args.lst:
        convert_lst(
            args.lst,
            not args.no_jacobian,
            args.qss_format_input,
            args.qss_symbolic_jacobian,
        )
    elif args.lst_qss:
        convert_lst_qss(
            args.lst_qss,
            not args.no_jacobian,
        )
    end = time.time()
    print(f"CEPTR run time: {end-start:.2f} s")


if __name__ == "__main__":
    main()
