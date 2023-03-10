import cantera as ct
import numpy as np
import sys
sys.path.append('../../utils')
import fileio as io
import os

gas1 = ct.Solution('mechanism/hychem.cti')
air = "O2:0.21,N2:0.79"
gas1.TPX = 1111, 1.0*ct.one_atm, 'O2: 0.207818, POSF11498: 0.0103909, N2: 0.781791'
gas1.transport_model='Mix'
#r = ct.IdealGasConstPressureReactor(gas1)
r = ct.IdealGasReactor(gas1)
sim = ct.ReactorNet([r])
time = 0.0
states = ct.SolutionArray(gas1, extra=['t'])
dt =  0.01
ndt = int(1000)

gas1()

timeTemp = np.zeros((ndt,5))
indexPOSF11498 = gas1.species_names.index('POSF11498')
indexO2 = gas1.species_names.index('O2')
indexCO2 = gas1.species_names.index('CO2')
for n in range(ndt):
    time += dt/ndt
    sim.advance(time)
    states.append(r.thermo.state, t=time*1e3)
    timeTemp[n,0]=time
    timeTemp[n,1]=r.T
    timeTemp[n,2]=r.Y[indexPOSF11498]
    timeTemp[n,3]=r.Y[indexO2]
    timeTemp[n,4]=r.Y[indexCO2]
print("Final Tem  = " , r.T)
print("Final Rho = ", r.density)
print("Final Y[POSF11498] = ", r.Y[indexPOSF11498])
print("Final Y[O2] = ", r.Y[indexO2])
print("Final Y[CO2] = ", r.Y[indexCO2])
print("h = ", r.thermo.h)

# Write result to file
os.makedirs('output',exist_ok=True)
io.writeMultiColFile('output/CanteraReaction_hychem.txt','',timeTemp)


