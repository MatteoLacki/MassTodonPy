import networkx as nx

class SFGraph(nx.Graph):
    def G_nodes(G):
        for G_name in G:
            if G.node[G_name]['type']=='G':
                yield G_name, G[G_name]


SFG = SFGraph()
SFG.add_node('M1', type='M')
SFG.add_node('I1', type='I')
SFG.add_node('I2', type='I')
SFG.add_node('I3', type='I')
SFG.add_node('I4', type='I')
SFG.add_node('G1', type='G')
SFG.add_node('G2', type='G')
SFG.add_edges_from([('M1','I1'),('M1','I2'),('M1','I3'),('M1','I4'),('G1','I1')])

for G_name, G in SFG.G_nodes():
    print G_name, G
