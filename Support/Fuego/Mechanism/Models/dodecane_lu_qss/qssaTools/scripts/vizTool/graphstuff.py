import networkx as nx
import matplotlib.pyplot as plt

G = nx.DiGraph()
G.add_edge('CH4','reac',weight=1)
G.add_edge('CO2','reac',weight=1)
G.add_edge('reac','H2O',weight=1)
G.add_edge('A','C',weight=1)
G.add_edge('E','C',weight=1)


# Specify the edges you want here
red_edges = [('A', 'C'), ('E', 'C')]
edge_colours = ['black' if not edge in red_edges else 'red'
                for edge in G.edges()]
black_edges = [edge for edge in G.edges() if edge not in red_edges]

nodeList = [node for node in G.nodes()]
nodeSize = [500 for node in nodeList]

# Need to create a layout when doing
# separate calls to draw nodes and edges
pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, node_size =nodeSize)
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos, edgelist=red_edges, edge_color='r', arrows=True)
nx.draw_networkx_edges(G, pos, edgelist=black_edges, arrows=True)
plt.show()
