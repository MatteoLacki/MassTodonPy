import igraph as ig
g = ig.Graph()
g.add_vertex( name='Dupa', big=1 )
g.add_vertex( name='Osiol', big=2 )
g.add_vertex( name='Chuj', big=2 )
g.add_vertex( name='Cipa', big=2 )
g.add_edge('Dupa','Osiol')
g.add_edge('Dupa','Cipa')

# print g.vs['name']
g.vs['gender'] = ['f','m','m','f']
color_dict = {'f': 'pink', 'm': 'blue'}
g.vs['color'] = [ color_dict[ge] for ge in g.vs['gender']]
print g

layout = g.layout("star")
g.vs['label'] = g.vs['name']
ig.plot(g, layout = layout, target='/Users/matteo/test1.pdf',  bbox = (300, 300), margin = 20)

g.delete_edges([('Dupa','Cipa')])
ig.plot(g, layout = layout, target='/Users/matteo/test2.pdf',  bbox = (300, 300), margin = 20)