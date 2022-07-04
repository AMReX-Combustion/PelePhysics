"""Convert cantera mechanism to C++ files."""
import argparse
import pathlib

import cantera as ct

import ceptr.converter as converter


def convert(fname, hformat, remove_1, remove_pow2, min_op_count):
    """Convert a mechanism file."""
    mechanism = ct.Solution(fname)
    conv = converter.Converter(
        mechanism, hformat, remove_1, remove_pow2, min_op_count
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
        "-rp2",
        "--remove_pow2",
        action="store_true",
        help="Remove pow(...,2) in printed expressions",
    )

    parser.add_argument(
        "-moc",
        "--min_op_count",
        type=int,
        metavar="",
        required=False,
        help="Min number of operation count per expression",
        default=0,
    )

    args = parser.parse_args()

    if args.fname:
        convert(
            args.fname,
            args.hformat,
            args.remove_1,
            args.remove_pow2,
            args.min_op_count,
        )
    elif args.lst:
        convert_lst(
            args.lst,
            args.hformat,
            args.remove_1,
            args.remove_pow2,
            args.min_op_count,
        )


if __name__ == "__main__":
    main()
