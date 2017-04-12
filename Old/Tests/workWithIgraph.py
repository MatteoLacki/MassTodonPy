import igraph as ig


V = [   {'name':('V1','ha'), 'mass':10, 'intensity':145},
        {'name':'V2', 'mass':30, 'intensity':1323},
        {'name':'V3', 'mass':2, 'intensity':1323}      ]

E = [   { 'source': ('V1','ha'), 'target': 'V2'},
        { 'source': ('V1','ha'), 'target': 'V3'}    ]

G_lists = ig.Graph.DictList(vertices=V, edges=E)

def plott( G, **kwds ):
	# elem_color 	= {'C':'grey', 'H':'red', 'N':'blue', 'O':'white', 'S':'pink' }
	visual_style= {}
	# visual_style['vertex_color']= [ elem_color[gender] for gender in G.vs['elem']]
	visual_style['vertex_label']= G.vs['name']
	visual_style.update(kwds)
	return ig.plot( G, **visual_style )

plott(G_lists)

a = V.__iter__()

G = ig.Graph.DictList( vertices = V.__iter__(), edges=E.__iter__() )
plott(G)

G.vs.degree()

