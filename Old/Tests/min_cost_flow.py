import networkx as nx
G = nx.DiGraph()
G.add_node('a', demand = -5)
G.add_node('d', demand = 5)
G.add_edge('a', 'b', weight = 3, capacity = 4)
G.add_edge('a', 'c', weight = 6, capacity = 10)
G.add_edge('b', 'd', weight = 1, capacity = 9)
G.add_edge('c', 'd', weight = 2, capacity = 5)
flowDict = nx.min_cost_flow(G)

flowDict

import  scipy as spy
import  pylab as pl
import  matplotlib.pyplot as plt

nx.draw(G,pos=nx.circular_layout(G))
plt.show()

Q = 8
qC= 2
qZ= 3
G = nx.DiGraph()
G.add_node('c', demand = -10)
G.add_node('z', demand = -15)
G.add_node('D', demand =  25)
G.add_node('cz_junction')

G.add_edge('c', 'D', weight = Q-1-qC)
G.add_edge('z', 'D', weight = Q-1-qZ)
G.add_edge('c', 'cz_junction')
G.add_edge('z', 'cz_junction')
G.add_edge('cz_junction', 'D', weight = Q-1-qC-qZ)

flowDict = nx.min_cost_flow(G)
flowDict

import  scipy as spy
import  pylab as pl
import  matplotlib.pyplot as plt

nx.draw(G,pos=nx.circular_layout(G), with_labels=True)
plt.show()

nx.min_cost_flow_cost(G)
