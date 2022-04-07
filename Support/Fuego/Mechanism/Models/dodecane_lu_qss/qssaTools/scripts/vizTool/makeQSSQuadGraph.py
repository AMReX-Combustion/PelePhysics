from FMC import FMC
import fuego
from QSSspecies import getListSpecies 
from coupling import *
import sys
sys.path.append('scripts/utils')
import myparser
import os


import networkx as nx
import matplotlib.pyplot as plt



# Parse input
inpt = myparser.parseInputFile()

# Set up filenames
mech_filename = inpt['mech']
therm_filename = inpt['therm']
tran_filename = inpt['tran']

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
species_non_qssa = getListSpecies(inpt['non_qssa_list'])
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

# Need to create a layout when doing
# separate calls to draw nodes and edges
pos = nx.kamada_kawai_layout(G)
##pos = nx.spring_layout(G)
##pos = nx.spiral_layout(G)

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
outputFolder = inpt['outputFolder']
plt.savefig(os.path.join(outputFolder,'quadGraphQSS.png'))
#plt.savefig(os.path.join(outputFolder,'quadGraphQSS.eps'))
#plt.show()



