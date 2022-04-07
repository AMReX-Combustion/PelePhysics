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
outputFolder = inpt['outputFolder'] 
plt.savefig(os.path.join(outputFolder,'directedGraphQSS.png'))
#plt.savefig(os.path.join(outputFolder,'directedGraphQSS.eps'))
#plt.show()



