# Given Epsilon in J/mol, calculate the Lennard-Jones well depth in Kelvin.

# Epsilon in J/mol
epsilon = 4954.04

def epsilon_to_lj_welldepth(epsilon_j_per_mol):
    """Convert Lennard-Jones epsilon from J/mol to K."""
    kb = 1.380649e-23  # Boltzmann constant in J/K
    Na = 6.02214076e23  # Avogadro's number in 1/mol
    epsilon_per_molecule = epsilon_j_per_mol / Na
    lj_welldepth = epsilon_per_molecule / kb
    print(f"Lennard-Jones well depth: {lj_welldepth:.3f} K") 

epsilon_to_lj_welldepth(epsilon)