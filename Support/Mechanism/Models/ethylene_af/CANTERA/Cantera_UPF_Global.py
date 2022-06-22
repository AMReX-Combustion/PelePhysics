###############################################################
#
# Cantera UPF : freely-propagating premixed flames
# Description : perform calculations of UPF for a range of eq. ratio
#               and compute global indices (integrated fuel, CO, NO consumption, ...)
# Author : A. Felden
# Last modified : 02/2018
#
###############################################################

#import :
from cantera import *
import numpy as np
import csv
import sys,os

#################################################################
# Prepare your run
#################################################################
#Mechanism used for the process
cas         = 'DET.cti' 
pref        = 'C2H4DET_T298_P1'
restore     = 'no'
init_flame  = 'RedHighP_HighP_save.xml' 

#General parameter values :
p          	= 101325.0                # pressure
tin        	= 298.0        # unburned gas temperature
phi_min 	= 0.8
phi_max 	= 1.5
npoints 	= 1

########
#STORAGE
########
# KSI
ksi = np.zeros(npoints,'d')
# YC eq
Yc_equil = np.zeros(npoints,'d')
#phiv
phi = np.zeros(npoints,'d')
#sl0phiv
sl_phi = np.zeros(npoints,'d')
#tadphiv
t_phi = np.zeros(npoints,'d')
#omega0phiv
#delta0phiv
dl0_phi = np.zeros(npoints,'d')

#########################################
# SPECS : profiles, wdot, consump, val max
# prof  : Fuel, oxy,  CO, CO2, OH, C2H2
# max   :       oxy,  CO, CO2, OH, C2H2
# wdot  : Fuel, oxy,  CO,          C2H2
# cons  : Fuel,       CO,          C2H2 
#########################################
#FUEL consump
Fuel_c 		= np.zeros(npoints,'d')
#OXY consump
O2m_phi 	= np.zeros(npoints,'d')
#CO consump
COm_phi 	= np.zeros(npoints,'d')
CO_c 		= np.zeros(npoints,'d')
#CO2 consump
CO2m_phi 	= np.zeros(npoints,'d')
#OH consump
OHm_phi 	= np.zeros(npoints,'d')


#################
#Run parameters:

#Initial grids, chosen to be 15cm long : 
initial_grid = np.array([0.0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1], 'd')/0.1 #, 0.11, 0.115, 0.119, 0.12],'d')/0.1 # m

#Tolerance properties
tol_ss    = [1.0e-8, 1.0e-12]        # [rtol atol] for steady-state problem
tol_ts    = [1.0e-8, 1.0e-12]        # [rtol atol] for time stepping

loglevel  = 1                       # amount of diagnostic output
				    
refine_grid = True                  # True to enable refinement

#################
#Stoechiometry :

fuel_species = 'C2H4'
stoich_O2 = 3.0
air_N2_O2_ratio = 3.76

gas  = Solution(cas)
m=gas.n_species

ifuel = gas.species_index(fuel_species)
ioh = gas.species_index('OH')
ih2 = gas.species_index('H2')
io2 = gas.species_index('O2')
in2 = gas.species_index('N2')
ico = gas.species_index('CO')
ico2 = gas.species_index('CO2')
ih2o = gas.species_index('H2O')

#################################################################
#BILGER PRE PROC :
#################################################################
#Mixt frac calculation based on Bilger
#molecular weight of CHON
C_W = 0.012
H_W = 0.001
O_W = 0.016
N_W = 0.014
#Bilger coefs pour formule
Bilger_coefs 		= np.zeros(4,'d')
Bilger_coefs[0]		= 2.0 / C_W
Bilger_coefs[1]		= 1.0 / ( 2.0 * H_W )
Bilger_coefs[2]		= -1.0 / O_W
Bilger_coefs[3]		= 0.0
#Storage nb d atomes par especes
Nb_C			= np.zeros(m,'d')
Nb_H			= np.zeros(m,'d')
Nb_O			= np.zeros(m,'d')
#Compo du fuel et oxi tank
Ftank			= np.zeros(m,'d')
Ftank[ifuel]		= 1.0
Oxtank			= np.zeros(m,'d')
Oxtank[io2]		= 0.233
Oxtank[in2]		= 0.767
#Normalisation
BetaF 			= 0.0
BetaOx 			= 0.0
for i, spec in enumerate(gas.species_names):
        Nb_C[i] = gas.n_atoms(spec,'C')  
        Nb_H[i] = gas.n_atoms(spec,'H')  
        Nb_O[i] = gas.n_atoms(spec,'O')  
        Weightedatom = Bilger_coefs[0] * C_W * Nb_C[i] + Bilger_coefs[1] * H_W * Nb_H[i] + Bilger_coefs[2] * O_W * Nb_O[i]
        BetaF = BetaF + 1000.0 * Weightedatom * Ftank[i] / gas.molecular_weights[i]
        BetaOx = BetaOx + 1000.0 * Weightedatom * Oxtank[i] / gas.molecular_weights[i]

if (restore == 'yes'): 
#################################################################
# RESTORE PREVIOUS FLAME
#################################################################
    phi[0]   = phi_min 
    x        = np.zeros(m,'d')
    x[ifuel] = phi[0]
    x[io2]   = stoich_O2
    x[in2]   = stoich_O2*air_N2_O2_ratio
    
    gas.TPX = tin, p, x
    
    f = FreeFlame(gas, initial_grid)
    f.restore(init_flame, 'phi_'+str(phi[0]))
else:
#################################################################
#FIRST SIMULATION :
#################################################################
#Stoechiometry :
    phi[0]   = phi_min 
    x        = np.zeros(m,'d')
    x[ifuel] = phi[0]
    x[io2]   = stoich_O2
    x[in2]   = stoich_O2*air_N2_O2_ratio

    print('  ')
    print(' ### First flame at Phi : '+str(phi[0]))
    
    gas.TPX = tin, p, x

#Create the free laminar premixed flame
    f = FreeFlame(gas, initial_grid)
    
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    
    f.inlet.X = x
    f.inlet.T = tin
    
    f.transport_model = 'Mix'

#First flame:

      #No energy for starters
    f.energy_enabled = False

      #Refinement criteria
    f.set_refine_criteria(ratio = 5.0, slope = 1, curve = 1)

      #Max number of times the Jacobian will be used before it must be re-evaluated
    f.set_max_jac_age(30, 30)

      #Set time steps whenever Newton convergence fails
    f.set_time_step(1.0e-6, [5, 10, 20]) #s

      #Calculation
    f.solve(loglevel, refine_grid)

#################
#Second flame:

	#Energy equation enabled
    f.energy_enabled = True

      #Refinement criteria when energy equation is enabled
    f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)

      #Calculation and save of the results
    f.solve(loglevel, refine_grid)

      #Refinement criteria when energy equation is enabled
    f.set_refine_criteria(ratio = 3.0, slope = 0.2, curve = 0.2)

      #Calculation and save of the results
    f.solve(loglevel, refine_grid)

#################
#Third flame and so on ...:

	#Refinement criteria should be changed ...
    f.set_refine_criteria(ratio = 2.0, slope = 0.04, curve = 0.04)

    f.solve(loglevel, refine_grid)

	#Refinement criteria should be changed ...
    f.set_refine_criteria(ratio = 2.0, slope = 0.02, curve = 0.02)

    f.solve(loglevel, refine_grid)

    sl_phi[0] = f.u[0]
    t_phi[0]  = f.T[f.flame.n_points-1]
    print( "Sl ? ",  sl_phi[0])
    print( "Tad ? ", t_phi[0])

	#...and saving of the results
    #f.save(str(pref)+'_save.xml','phi_'+str(phi[0]))
    f.write_csv(str(pref)+'_Phi'+str(phi[0])+'_Y.csv', species='Y', quiet=False)

