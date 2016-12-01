import igraph as ig

# g = ig.Graph( directed = True )
g = ig.Graph()
g.add_vertices( 5 )
g.add_edges( [(0,1), (0,2),(0,3), (0,4), (1,0)] )
g.vs['type'] = ['f','p','p','p','p']


# ig.to_directed(g, mode = 'mutual')
g = g.as_directed('mutual')


# h = ig.Graph( directed = True )
# h.add_vertices( 3 )
# h.add_edges( [(0,1), (0,2)] )

# j = g+h
# dupa = j.components(mode='strong')

# # print dupa
# for d in dupa:
# 	print d


h = ig.Graph( directed=True )
h.add_vertices( 6 )
h.add_edges( [(0,1),(1,3),(3,1),(1,4),(4,1),(0,2),(2,5),(5,2)] )

# print h.components(mode='strong')
# print 

h.es['ban'] = ['a','b','c','d','e','f','g','h']
h.vs['name']= range(6,0,-1) # There are 6 vertices.
# print h.vs['invNo']

# if h[0,1]:
# 	print 'hello'
# 	print h
# else:
# 	print 'dupa'
# print h.es[0]
# print h.es[[(0,1)]]
# print
# print h.vs['invNo']
# h.delete_edges([(0,1)])
# print h
# print h.es[0]
# h.delete_edges(0)
# print h.es[0]

print h
# h.delete_vertices()
# print h
layout = h.layout("kk")
ig.plot(h, layout = layout)