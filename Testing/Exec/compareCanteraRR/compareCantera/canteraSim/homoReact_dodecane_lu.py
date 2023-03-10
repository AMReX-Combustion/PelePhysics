import cantera as ct
import numpy as np
import sys
sys.path.append('../../utils')
import fileio as io
import os

gas1 = ct.Solution('mechanism/dodecane_lu.cti')
gas1.TPX = 600.0, 5.0*ct.one_atm, 'O2: 0.20153, NC12H26: 0.040307, N2: 0.75816'
gas1.transport_model='Mix'
#r = ct.IdealGasConstPressureReactor(gas1)
r = ct.IdealGasReactor(gas1)
sim = ct.ReactorNet([r])
time = 0.0
states = ct.SolutionArray(gas1, extra=['t'])
dt =  3.0
ndt = int(1000)

gas1()

timeTemp = np.zeros((ndt,5))
indexNC12H26 = gas1.species_names.index('NC12H26')
indexO2 = gas1.species_names.index('O2')
indexCO2 = gas1.species_names.index('CO2')
for n in range(ndt):
    time += dt/ndt
    sim.advance(time)
    states.append(r.thermo.state, t=time*1e3)
    timeTemp[n,0]=time
    timeTemp[n,1]=r.T
    timeTemp[n,2]=r.Y[indexNC12H26]
    timeTemp[n,3]=r.Y[indexO2]
    timeTemp[n,4]=r.Y[indexCO2]
print("Final Tem  = " , r.T)
print("Final Rho = ", r.density)
print("Final Y[NC12H26] = ", r.Y[indexNC12H26])
print("Final Y[O2] = ", r.Y[indexO2])
print("Final Y[CO2] = ", r.Y[indexCO2])
print("h = ", r.thermo.h)

# Write result to file
os.makedirs('output',exist_ok=True)
io.writeMultiColFile('output/CanteraReaction_dodecane_lu.txt','',timeTemp)


