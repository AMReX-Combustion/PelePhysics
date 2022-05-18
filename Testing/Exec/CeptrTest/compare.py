"""Compare Cantera and PelePhysics results."""
import argparse

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

plt.rc("text", usetex=True)
cmap_med = [
    "#F15A60",
    "#7AC36A",
    "#5A9BD4",
    "#FAA75B",
    "#9E67AB",
    "#CE7058",
    "#D77FB4",
    "#737373",
]
cmap = [
    "#EE2E2F",
    "#008C48",
    "#185AA9",
    "#F47D23",
    "#662C91",
    "#A21D21",
    "#B43894",
    "#010202",
]
dashseq = [
    (None, None),
    [10, 5],
    [10, 4, 3, 4],
    [3, 3],
    [10, 4, 3, 4, 3, 4],
    [3, 3],
    [3, 3],
]
markertype = ["s", "d", "o", "p", "h"]

attrs = {"temperature": {"label": r"$T~[\mathrm{K}]$"}}


def main():
    """Compare results."""
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-c", "--cname", help="Cantera results file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--pname", help="PelePhysics results file", type=str, required=True
    )
    args = parser.parse_args()

    cdf = pd.read_csv(args.cname)
    pdf = pd.read_csv(args.pname)

    fields = ["temperature"]
    for field in fields:
        plt.figure(field)
        p = plt.plot(cdf.time, cdf[field], lw=2, color=cmap[0], label="Cantera")
        p[0].set_dashes(dashseq[0])
        p = plt.plot(pdf.time, pdf[field], lw=2, color=cmap[1], label="PelePhysics")
        p[0].set_dashes(dashseq[1])

    fname = "plots.pdf"
    with PdfPages(fname) as pdf:
        for field in fields:
            plt.figure(field)
            ax = plt.gca()
            plt.xlabel(r"$t~[\mathrm{s}]$", fontsize=22, fontweight="bold")
            plt.ylabel(attrs[field]["label"], fontsize=22, fontweight="bold")
            plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
            plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
            ax.legend(loc="best")
            plt.tight_layout()
            pdf.savefig(dpi=300)


if __name__ == "__main__":
    main()
