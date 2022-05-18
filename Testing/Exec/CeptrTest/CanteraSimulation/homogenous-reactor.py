"""Run a homegenous reactor in Cantera."""
import argparse

import cantera as ct
import numpy as np


def main():
    """Run the reactor."""
    parser = argparse.ArgumentParser(description="Cantera homogenous reactor")
    parser.add_argument("-f", "--fname", help="Mechanism file", type=str, required=True)
    args = parser.parse_args()

    mechanism = ct.Solution(args.fname)
    mechanism.TPX = 600.0, 5 * ct.one_atm, "H2:0.0,O2:0.7,NC12H26:0.3"
    mechanism.transport_model = "Mix"
    # r = ct.IdealGasConstPressureReactor(mechanism)
    r = ct.IdealGasReactor(mechanism)
    sim = ct.ReactorNet([r])
    time = 0.0
    states = ct.SolutionArray(mechanism, extra=["t"])
    dt = 10
    ndt = int(1000)

    mechanism()

    time_temp = np.zeros((ndt, 2))
    for n in range(ndt):
        time += dt / ndt
        sim.advance(time)
        states.append(r.thermo.state, t=time * 1e3)
        time_temp[n, 0] = time
        time_temp[n, 1] = r.T

    # Write result to file
    oname = "cantera_results.txt"
    header = ["time", "temperature"]
    np.savetxt(oname, time_temp, delimiter=",", header=",".join(header), comments="")


if __name__ == "__main__":
    main()
