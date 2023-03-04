"""
Constant-pressure, adiabatic kinetics simulation.

Requires: cantera >= 2.5.0, matplotlib >= 2.0
Compares the ignition delay time from the original mechanism
to the mechanism with the PAH module added
"""

import sys

import numpy as np
import cantera as ct
import matplotlib.pyplot as plt


def compareIDT(origfile, newfile, numatm, fuel):
    files = [origfile, newfile]
    pres = numatm * ct.one_atm
    title = fuel + "_"+str(int(numatm)) + "atm"
    print("Comparing IDT at " + str(numatm) + " atm")
    if (numatm == 1):
        Tlo = 1000.
        Thi = 1600.
    else:
        Tlo = 400.
        Thi = 1250.
    Tvals = np.arange(Tlo, Thi, 50.)
    phivals = [0.5, 1., 1.3]

    def getIDT(time, T, threshold):
        # Get ignition delay time
        if (np.max(T) < threshold):
            return 0.
        else:
            above = T >= threshold
            IDT = np.min(time[above])
            return IDT
    def getIDTFirst(time, X):
        threshold = np.max(X)
        above = X >= threshold
        IDT = np.min(time[above])
        return IDT
    def solveIDT(filename, T0, pres, equiv, fuel):
        delta_T_max = 20.
        dt_max = 1.e-5
        t_end = 200.*dt_max
        gas = ct.Solution(filename)
        gas.TP = T0, pres
        gas.set_equivalence_ratio(equiv, fuel=fuel, oxidizer="O2:0.233, N2:0.767",basis="mass")
        r = ct.IdealGasConstPressureReactor(gas)
        sim = ct.ReactorNet([r])
        sim.verbose = False
        # limit advance when temperature difference is exceeded
        r.set_advance_limit('temperature', delta_T_max)
        states = ct.SolutionArray(gas, extra=['t'])
        while sim.time < t_end:
            sim.advance(sim.time + dt_max)
            states.append(r.thermo.state, t=sim.time*1e3)
        time = np.array(states.t)
        T = np.array(states.T)
        heat = np.array(states.heat_release_rate)
        h2o2 = np.array(states.Y[:, gas.species_index("H2O2")])
        idt0 = getIDT(time, T, 1600.)
        idt1 = getIDTFirst(time, h2o2)
        return [idt0, idt1]

    fig, ax = plt.subplots()
    for phi in phivals:
        print("Comparing at equivalence ratio " + str(phi))
        aTval = []
        idtval0 = []
        idtval1 = []
        for T in Tvals:
            idt = np.zeros(len(files))
            idt1 = np.zeros(len(files))
            for cf in range(len(files)):
                curfile = files[cf]
                [idt0, idt2] = solveIDT(curfile, T, pres, phi, fuel)
                idt[cf] = idt0
                idt1[cf] = idt2
            if (idt[0] > 0.):
                aTval.append(T)
                idtval0.append(idt[0])
                idtval1.append(idt[1])
                pdiff = (idt[0] - idt[1]) / idt[0]
                print("T {:1.1f} phi {:1.2f} IDT values {:1.4f} {:1.4f} pdiff {:1.5f}".format(T, phi, idt[0], idt[1], pdiff))
        atitle = title + "_phi_" + str(phi)
        plt.figure()
        plt.title(atitle)
        atitle = atitle.replace(".", "")
        plt.scatter(aTval, idtval0, color='black', label="Orig")
        plt.scatter(aTval, idtval1, color='red', label="PAH mod")
        plt.xlabel("T (K)")
        plt.ylabel("Ign. Delay (ms)")
        plt.legend()
        plt.savefig(atitle + ".png")


