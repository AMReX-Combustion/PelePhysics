import cantera as ct
import numpy as np
import sys
sys.path.append('../../utils')
import fileio as io
import os

gas1 = ct.Solution('mechanism/dodecane_lu.cti')
gas1.TPX = 600.0, 5*ct.one_atm, 'H2:0.0,O2:0.7,NC12H26:0.3'
gas1.transport_model='Mix'
#r = ct.IdealGasConstPressureReactor(gas1)
r = ct.IdealGasReactor(gas1)
sim = ct.ReactorNet([r])
time = 0.0
states = ct.SolutionArray(gas1, extra=['t'])
dt =  2
ndt = int(10000)

gas1()

timeTemp = np.zeros((ndt,8))
#print('%10s %10s %10s %14s' % ('t [s]','T [K]','P [Pa]','u [J/kg]'))
indexNC12H26 = gas1.species_names.index('NC12H26')
indexH2 = gas1.species_names.index('H2')
indexO2 = gas1.species_names.index('O2')
indexCO2 = gas1.species_names.index('CO2')
indexCO = gas1.species_names.index('CO')
for n in range(ndt):
    time += dt/ndt
    sim.advance(time)
    states.append(r.thermo.state, t=time*1e3)
    #print('%10.3e %10.3f %10.3f %14.6e' % (sim.time, r.T,
    #                                       r.thermo.P, r.thermo.u))
    timeTemp[n,0]=time
    timeTemp[n,1]=r.T
    timeTemp[n,2]=r.Y[indexNC12H26]
    timeTemp[n,3]=r.Y[indexO2]
    timeTemp[n,4]=r.Y[indexCO2]
    timeTemp[n,5]=r.Y[indexCO]
    timeTemp[n,6]=r.Y[indexH2]
    timeTemp[n,7]=r.thermo.h
print("Final Tem  = " , r.T)
print("Final Rho = ", r.density)
print("Final mass = ", r.mass)
print("Final Y[NC12H26] = ", r.Y[indexNC12H26])
print("Final Y[H2] = ", r.Y[indexH2])
print("Final Y[O2] = ", r.Y[indexO2])
print("h = ", r.thermo.h)

# Write result to file
os.makedirs('output',exist_ok=True)
io.writeMultiColFile('output/CanteraReaction_dodecane.txt','',timeTemp)


