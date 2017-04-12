import networkx as nx

G = nx.Graph()



# Any hashable object can become a node
G.add_node('Amphetamine')
G.add_node(('Amphetamine',100))

G.nodes()
G.node['Amphetamine']

G.add_nodes_from('ABC')

G.node['A']['bla']=10
G.nodes()
G.node['Amphetamine']
['blabla'] = 12
|H = nx.path_graph(10)
G.add_nodes_from(H)

for g in G.nodes_iter():
    print g

def stupidIter(N,W):
    for i in xrange(W,N+W):
        yield i

G.add_nodes_from(stupidIter(100,40))

G.neighbors('Amphetamine')
G.nodes()
G.edges()

for e in G.edges_iter():
    print e

for g in G.nodes_iter():
    print g

G.add_edge(1,3)
G[1][3]['color'] = 'blue'
G[3][1]

############################################################
FG = nx.Graph()
FG.add_weighted_edges_from([(1,2,0.125),(1,3,0.75),(2,4,1.2),(3,4,0.375)])

FG.nodes()
FG.edges()

for e in FG.edges_iter():
    print e

for n, nbrs in FG.adjacency_iter():
    for nbr,eattr in nbrs.items():
        data = eattr['weight']
        if data<0.5: print('(%d, %d, %.3f)' % (n,nbr,data))
FG.graph['day'] = 'Monday'
FG.add_node('Hurray', time='5pm')
FG.nodes(data=True)
############################################################
