import networkx as nx
G = nx.DiGraph()
G.add_edge('x','a', capacity=3.0)
G.add_edge('x','b', capacity=1.0)
G.add_edge('a','c', capacity=3.0)
G.add_edge('b','c', capacity=5.0)
G.add_edge('b','d', capacity=4.0)
G.add_edge('d','e', capacity=2.0)
G.add_edge('c','y', capacity=2.0)
G.add_edge('e','y', capacity=3.0)

import  scipy as spy
import  pylab as pl
import  matplotlib.pyplot as plt

nx.draw(G,pos=nx.circular_layout(G))
plt.show()

flow_value = nx.maximum_flow_value(G, 'x', 'y')
flow_value
