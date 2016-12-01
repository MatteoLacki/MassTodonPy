import igraph as ig

g = ig.Graph()
g.add_vertices(10)
g.add_edges([(0,1), (1,2),(2,3), (4,5)])
print g.components()
print g


h = ig.Graph()
h.add_vertices(3)
h.add_edges([(0,2)])
print h
z = g+h

print z.components()