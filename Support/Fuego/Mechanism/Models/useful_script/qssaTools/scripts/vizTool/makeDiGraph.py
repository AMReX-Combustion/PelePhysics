from FMC import FMC
import argparse
import fuego
from QSSspecies import getListSpecies 
from coupling import *
import sys
import os
try:
    import networkx as nx
except ImportError:
    sys.exit('Need networkx python module')
try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit('Need matplotlib python module')

# CLI
parser = argparse.ArgumentParser(description='Plot QSS species dependency graph')
parser.add_argument('-m', '--mechanism', type=str, metavar='', required=False, help='QSS mechanism without coupling treatment', default='output/qssa_nostar.inp')
parser.add_argument('-th', '--thermoFile', type=str, metavar='', required=False, help='Thermodynamic file', default='output/therm_nostar.dat')
parser.add_argument('-tr', '--tranFile', type=str, metavar='', required=False, help='Transport file', default='output/tran_nostar.dat')
parser.add_argument('-nqss', '--nonQSSSpecies', type=str, metavar='', required=False, help='Filename with non QSS species names', default='output/non_qssa_list_nostar.txt')
parser.add_argument('-o', '--outputFolder', type=str, metavar='', required=False, help='Where to store by product of the preprocessing',default='output')
group = parser.add_mutually_exclusive_group()
group.add_argument('-q','--quiet', action='store_true', help='Execute without plotting on screen')
group.add_argument('-v','--verbose', action='store_true', help='Plot on screen')
args = parser.parse_args()

# Set up filenames
mech_filename = args.mechanism
therm_filename = args.thermoFile
tran_filename = args.tranFile
non_qssa_species = args.nonQSSSpecies
outputFolder = args.outputFolder

# Load skeletal mechanism
app = FMC()
app.inventory.mechanism = mech_filename
app.inventory.thermo = therm_filename
app.inventory.trans = tran_filename
mechanism = fuego.serialization.mechanism()
mechanism.externalThermoDatabase(app.inventory.thermo)
mechanism.externalTransDatabase(app.inventory.trans)
mechanism = fuego.serialization.load(filename=app.inventory.mechanism, format='chemkin', mechanism=mechanism)
reactionIndex = mechanism._sort_reactions()

# Get list of intended QSS species
species_non_qssa = getListSpecies(non_qssa_species)
species_all = getListSpecies(app.inventory.mechanism)
species_qssa = list(set(species_all) - set(species_non_qssa))


# Figure out links between QSS species
links = getDiGraphLinks(mechanism,species_qssa)


# Make Graph
G = nx.DiGraph()
G.add_edges_from(links)

# Specify the edges you want here
nodeSize = []
labels = {}  
# define color map. user_node = red, book_nodes = green
nodeColor = []        
for node in G.nodes():
     nodeSize.append(0)
     labels[node] = node
     nodeColor.append('black')

black_edges = [edge for edge in G.edges()]

# Need to create a layout when doing
# separate calls to draw nodes and edges
#pos = nx.kamada_kawai_layout(G)
pos = nx.spring_layout(G)
##pos = nx.spiral_layout(G)
nx.draw_networkx_nodes(G, pos, node_size =nodeSize, node_color=nodeColor)
nx.draw_networkx_labels(G, pos,labels, font_size=14,font_color='k')
nx.draw_networkx_edges(G, pos, edgelist=black_edges, arrows=True)
plt.tight_layout()
plt.savefig(os.path.join(outputFolder,'directedGraphQSS.png'))


if args.verbose:
    plt.show()



