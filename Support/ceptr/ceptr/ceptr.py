"""Convert cantera mechanism to C++ files."""
import argparse
import pathlib

import cantera as ct

import ceptr.converter as converter


def convert(fname, hformat):
    """Convert a mechanism file."""
    mechanism = ct.Solution(fname)
    conv = converter.Converter(mechanism, hformat)
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
    parser = argparse.ArgumentParser(description="Mechanism converter")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fname", help="Mechanism file", type=str)
    group.add_argument(
        "-l", "--lst", help="Mechanism directory file list", type=str
    )

    parser.add_argument(
        "--hformat",
        help="sytle format for .H file output",
        type=str,
        choices=["readable", "gpu"],
        default="readable",
        required=False,
    )

    args = parser.parse_args()

    if args.fname:
        convert(args.fname, args.hformat)
    elif args.lst:
        convert_lst(args.lst, args.hformat)


if __name__ == "__main__":
    main()
