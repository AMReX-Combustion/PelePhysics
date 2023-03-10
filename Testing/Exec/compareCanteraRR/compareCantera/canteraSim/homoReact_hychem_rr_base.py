import cantera as ct
import numpy as np
import sys
import os

def printNonZeros(A, name):
    print(f"{name} :")
    for i in range(len(A)):
        if abs(A[i])>1e-300:
            #print(f"\t{i} : {A[i]}")
            print("\t"+str(i)+" : %.6g"%A[i])

gas1 = ct.Solution('mechanism/hychem.cti')
gas1.TPX = 900.0, 1.4*ct.one_atm, 'O2: 0.20153, POSF11498: 0.040307, N2: 0.75816'
gas1.transport_model='Mix'
#r = ct.IdealGasConstPressureReactor(gas1)
r = ct.IdealGasReactor(gas1)
sim = ct.ReactorNet([r])
time = 0.0
states = ct.SolutionArray(gas1, extra=['t'])
dt =  3.0
ndt = int(1000)

#breakpoint()

printNonZeros(gas1.concentrations*np.float64(1e-3), "sc")
printNonZeros(gas1.net_production_rates*np.float64(1e-9), "wdot")
printNonZeros(gas1.forward_rates_of_progress*np.float64(1e-9), "qf")
printNonZeros(gas1.reverse_rates_of_progress*np.float64(1e-9), "qr")

