###############################################################
# ADIABATIC_FLAME - A freely-propagating, premixed flat flame 
#                   of fuel in AIR (O2 + N2) 
###############################################################

#import :

from cantera import *
import numpy as np
import csv

#################################################################
# Prepare your run
#################################################################
# Parameter values :

# Mixture
mechanism        = 'drm19'
fuel_species     = 'CH4'
fuel_species_dil = 'H2'
restore          = 'yes'
	
# General
p            = 101325           # pressure [Pa]
tin          = 300.0            # unburned gas temperature [K]
phi          = 1.0              # Eq. ratio [-] 
dilution     = 0.0              # Dil by H2

# Refined grid at inlet and outlet, 6 points in x-direction :
domain_size  = 0.04             # Domain size [m] 
initial_grid = domain_size*np.array([0.0, 0.001, 0.01, 0.02, 0.029, 0.03],'d')/0.03 # m

#Set tolerance properties
tol_ss    = [1.0e-8, 1.0e-9]        # [rtol atol] for steady-state problem
tol_ts    = [1.0e-8, 1.0e-9]        # [rtol atol] for time stepping

loglevel  = 1                       # amount of diagnostic output
				    
refine_grid = True                  # True to enable refinement

# Print information
specformat = 'mole'
#POSSIBLY MODIFY THIS IF RESTORE IN XML WITH KNOWN NAME
if (dilution > 0.0):
    label_xml = 'FFdil-'+mechanism 
    label     = fuel_species+'_PHI'+str(phi)+'_T'+str(tin)+'_P'+str(p)
else:
    label_xml = 'FF-'+mechanism 
    label     = fuel_species+'_PHI'+str(phi)+'_T'+str(tin)+'_P'+str(p)+'_'+fuel_species_dil+'dil'+str(dilution)

#################
#Stoechiometry :
gas                = Solution(mechanism + '.cti','gas')
m                  = gas.n_species

stoich_fuel_O2     = gas.n_atoms(fuel_species,'C') + 0.25*gas.n_atoms(fuel_species,'H')
stoich_dil_O2      = 0.25*gas.n_atoms(fuel_species_dil,'H') 
stoich_hyb         = (1.0 - dilution) * stoich_fuel_O2 + dilution * stoich_dil_O2
air_N2_O2_ratio    = 3.76

# specs indx
ifuel = gas.species_index(fuel_species)
idil  = gas.species_index(fuel_species_dil)
io2   = gas.species_index('O2')
in2   = gas.species_index('N2')

# Mixture composition
x = np.zeros(m,'d')
x[ifuel]  = phi * (1.0 - dilution) 
x[idil]   = phi * dilution
x[io2]    = stoich_hyb
x[in2]    = stoich_hyb*air_N2_O2_ratio


#################
#Assembling objects :
if (restore == 'yes'): 

    #Set gas state to that of the unburned gas
    gas.TPX = tin, p, x
    
    #Create the free laminar premixed flame
    f = FreeFlame(gas, initial_grid)

    f.restore(label_xml+'.xml',label)

else:
    #Set gas state to that of the unburned gas
    gas.TPX = tin, p, x
    
    #Create the free laminar premixed flame
    f = FreeFlame(gas, initial_grid)
    
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    
    f.inlet.X = x
    f.inlet.T = tin
    
    f.transport_model = 'Mix'
    
    #################################################################
    # Program starts here
    #################################################################
    #First flame:
    
    #No energy for starters
    f.energy_enabled = False
    
    #Refinement criteria
    f.set_refine_criteria(ratio = 7.0, slope = 1, curve = 1)
    
    #Max number of times the Jacobian will be used before it must be re-evaluated
    f.set_max_jac_age(10, 10)
    
    #Set time steps whenever Newton convergence fails
    f.set_time_step(1.e-06, [1, 2, 5, 10]) #s
    f.max_time_step_count = 3000
    f.max_grid_points = 1000
    
    #Calculation
    f.solve(loglevel, refine_grid)
    
    #################
    #Second flame:
    
    #Energy equation enabled
    f.energy_enabled = True
    
    #Refinement criteria when energy equation is enabled
    f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)
    
    #Calculation 
    f.solve(loglevel, refine_grid)
    
    #################
    #Third flame and so on ...:
    f.set_refine_criteria(ratio = 3.0, slope = 0.1, curve = 0.1)
    
    f.solve(loglevel, refine_grid)
    
    ##################
    ##Fourth flame and so on ...:
    f.set_refine_criteria(ratio = 2.0, slope = 0.05, curve = 0.05, prune = 0.01)
    
    f.solve(loglevel, refine_grid)
    
    ##################
    ##Fifth flame and so on ...
    f.set_refine_criteria(ratio = 2.0, slope = 0.02, curve = 0.02, prune = 0.01)
    
    f.solve(loglevel, refine_grid)
    
    print('mixture averaged flamespeed = ',f.u[0])

    #Salve in xml format for later restart
    f.save(label_xml+'.xml',label)

#################################################################
# Save your results if needed
#################################################################
#Write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv(label_xml +'-'+ label+'.csv', quiet=True)

# PeleLM PMF file: MKS and mass fractions !
if (specformat == 'mole'):
    nz = f.flame.n_points
    csv_file = str("pmf"+label_xml +'-'+ label+"X.txt")
    with open(csv_file, 'w') as outfile:
        writer = csv.writer(outfile, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['VARIABLES = "X" "temp" "u" "rho"']+['"X_'+x+'"' for x in gas.species_names])
        writer.writerow(['ZONE I=' + str(nz) + ' FORMAT=POINT' + ' SPECFORMAT=MOLE'])
        for kk in range(len(f.grid)):
            f.set_gas_state(kk)
            writer.writerow([f.grid[kk],gas.T, f.u[kk], gas.density]+list(gas.X))
if (specformat == 'mass'):
    nz = f.flame.n_points
    csv_file = str("pmf"+label_xml +'-'+ label+"Y.txt")
    with open(csv_file, 'w') as outfile:
        writer = csv.writer(outfile, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['VARIABLES = "X" "temp" "u" "rho"']+['"Y_'+x+'"' for x in gas.species_names])
        writer.writerow(['ZONE I=' + str(nz) + ' FORMAT=POINT' + ' SPECFORMAT=MASS'])
        for kk in range(len(f.grid)):
            f.set_gas_state(kk)
            writer.writerow([f.grid[kk],gas.T, f.u[kk], gas.density]+list(gas.Y))
