"""Convert cantera mechanism to C++ files."""

import argparse
import pathlib
import time
from itertools import repeat
from multiprocessing import Pool, cpu_count

import cantera as ct

import ceptr.converter as converter


def parse_lst_file(lst):
    """Return mechanism paths give a file containing a list of mechanism files."""
    lpath = pathlib.Path(lst)
    fnames = []
    with open(lst, "r") as f:
        for line in f:
            if not line.startswith("#"):
                fnames.append(line)
    return [lpath.parents[0] / fn.strip() for fn in fnames]


def parse_qss_lst_file(lst):
    """Return mechanism paths give a file containing a list of qss mechanism files."""
    lpath = pathlib.Path(lst)
    fnames = []
    with open(lst, "r") as f:
        for line in f:
            if not line.startswith("#"):
                fnames.append(line)

    mechnames = [lpath.parents[0] / fn.split()[0].strip() for fn in fnames]
    qss_format_inputs = [lpath.parents[0] / fn.split()[1].strip() for fn in fnames]
    return mechnames, qss_format_inputs


def convert(
    fname,
    jacobian,
    qss_format_input,
    qss_symbolic_jac,
):
    """Convert a mechanism file."""
    print(f"""Converting file {fname}""")
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
    ncpu,
):
    """Convert mechanisms from a file containing a list of directories."""
    mechnames = parse_lst_file(lst)
    print(f"Using {ncpu} processes")
    with Pool(ncpu) as pool:
        pool.starmap(
            convert,
            zip(
                mechnames,
                repeat(jacobian),
                repeat(qss_format_input),
                repeat(qss_symbolic_jac),
            ),
        )


def convert_lst_qss(
    lst,
    jacobian,
    ncpu,
):
    """Convert QSS mechanisms from a file of directories and format input."""
    mechnames, qss_format_inputs = parse_qss_lst_file(lst)
    print(f"Using {ncpu} processes")
    with Pool(ncpu) as pool:
        pool.starmap(
            convert,
            zip(mechnames, repeat(jacobian), qss_format_inputs, repeat(True)),
        )


def main():
    """Convert cantera mechanisms to C++ files."""
    start = time.time()
    parser = argparse.ArgumentParser(
        description="Mechanism converter",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fname", help="Mechanism file", type=str)
    group.add_argument("-l", "--lst", help="Mechanism directory file list", type=str)
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

    parser.add_argument(
        "-n", "--ncpu", help="Number of processes to use", type=int, default=cpu_count()
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
            args.ncpu,
        )
    elif args.lst_qss:
        convert_lst_qss(
            args.lst_qss,
            not args.no_jacobian,
            args.ncpu,
        )
    end = time.time()
    print(f"CEPTR run time: {end-start:.2f} s")


if __name__ == "__main__":
    main()
