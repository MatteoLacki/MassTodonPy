import networkx as nx

H = nx.Graph()
H.add_edges_from([(1,2),(1,3)])
H.nodes()
H.remove_edge(1,2)
H.edges()
