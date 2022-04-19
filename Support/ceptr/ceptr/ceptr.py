"""Convert cantera mechanism to C++ files."""
import argparse

import cantera as ct

import ceptr.converter as converter


def main():
    """Convert cantera mechanism to C++ files."""
    parser = argparse.ArgumentParser(description="Mechanism converter")
    parser.add_argument(
        "-f", "--fname", help="Mechanism file", required=True, type=str
    )
    args = parser.parse_args()

    mechanism = ct.Solution(args.fname)
    conv = converter.Converter(mechanism)
    conv.writer()
    conv.formatter()


if __name__ == "__main__":
    main()
