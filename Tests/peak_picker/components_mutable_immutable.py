import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter



SFG.add_node('A', intensity=10., type='E', mz=0 )
SFG.add_node('B', intensity=10., type='E', mz=1 )
SFG.add_edge('A','B')

SFG.add_node('C', intensity=11., type='E', mz=2 )
SFG.add_node('D', intensity=11., type='E', mz=2 )
SFG.add_edge('C','D')


A, B = tuple(nx.connected_component_subgraphs(SFG))
A, B = tuple(nx.connected_components(SFG))

A.nodes(data=True)
A.node['A']['intensity']= 25
A.nodes(data=True)
SFG.nodes(data=True)
