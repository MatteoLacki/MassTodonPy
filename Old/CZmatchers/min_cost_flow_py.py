import networkx as nx

BFG = nx.DiGraph()
BFG.add_nodes_from(range(6))
BFG.add_edges_from([
                   (0, 2, {'capacity':4}),
                   (0, 4, {'capacity':5}),
                   (4, 3, {'cost':10.1}),
                   (2, 3, {'cost':5.10}),
                   (4, 5, {'cost':10.1}),
                   (3, 1, {'capacity':5}),
                   (5, 1, {'capacity':4})
                  ])

BFGC = BFG.copy()

print BFG.edges()

BFGC.node[0]['demand'] = -1000
BFGC.node[1]['demand'] = 1000
BFGC.add_edge(0, 1)
for v1, v2 in BFGC.edges():
   BFGC[v1][v2]['weight'] = -BFGC[v1][v2].get('cost', 0)

BFGC.edges(data=True)
print nx.min_cost_flow(BFGC)
