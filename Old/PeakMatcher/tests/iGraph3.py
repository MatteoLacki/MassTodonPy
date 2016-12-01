from igraph import *
# import igraph as ig

# g = ig.Graph( directed = True )
g = Graph()
g.add_vertices( 5 )
g.add_edges( [(0,1), (0,2),(0,3), (0,4), (1,0)] )
g.vs['type'] = ['f','p','p','p','p']

print g.vs['type']


# layout = g.layout("kk")
# plot(g, layout = layout)