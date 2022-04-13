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
parser = argparse.ArgumentParser(description='Highlight possible quadratic coupling')
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
links, quadraticReac, problematicSpecies = getQuadGraphLinks(mechanism,species_qssa)


# Make Graph
G = nx.DiGraph()
G.add_edges_from(links)

# Specify the edges you want here
labels = {}  
labels_normalSpecies = {}  
labels_problemSpecies = {}  
# define color map. user_node = red, book_nodes = green
nodeColor = []        
nodeSize = []        
node_list_normalSpecies=[]
node_list_problemSpecies=[]
node_list_normalReac=[]
node_list_problemReac=[]
for node in G.nodes():
    if node in quadraticReac:
        node_list_problemReac.append(node)
    elif node in problematicSpecies:
        node_list_problemSpecies.append(node)
        labels_problemSpecies[node] = node
    elif node.startswith('F') or node.startswith('R'):
        node_list_normalReac.append(node)
    else:
        node_list_normalSpecies.append(node)
        labels_normalSpecies[node] = node


black_edges = [edge for edge in G.edges()]

# Create a layout
pos = nx.kamada_kawai_layout(G)

# Draw normal species
nx.draw_networkx_nodes(G, pos, nodelist=node_list_normalSpecies, node_size=0, node_color='b')
label_options = {"ec": "k", "fc": 'white', "alpha": 1}
nx.draw_networkx_labels(G, pos, labels=labels_normalSpecies, nodeList=node_list_normalSpecies, font_size=14,font_color='k',bbox=label_options)

# Draw problem species
nx.draw_networkx_nodes(G, pos, nodelist=node_list_problemSpecies, node_size=0, node_color='b')
label_options = {"ec": "k", "fc": 'red', "alpha": 1}
nx.draw_networkx_labels(G, pos, labels=labels_problemSpecies, nodeList=node_list_problemSpecies, font_size=14,font_color='k',bbox=label_options)

# Draw normal reactions
nx.draw_networkx_nodes(G, pos, nodelist=node_list_normalReac, node_size=200, node_color='w', edgecolors='k')

# Draw problem reactions
nx.draw_networkx_nodes(G, pos, nodelist=node_list_problemReac, node_size=200, node_color='r')



nx.draw_networkx_edges(G, pos, edgelist=black_edges, arrows=False)
plt.tight_layout()
plt.savefig(os.path.join(outputFolder,'quadGraphQSS.png'))

if args.verbose:
    plt.show()

