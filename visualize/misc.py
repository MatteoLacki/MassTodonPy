import igraph as ig

def plott( G, **kwds ):
	elem_color 	= {'C':'grey', 'H':'red', 'N':'blue', 'O':'white', 'S':'pink' }
	visual_style= {}
	visual_style['vertex_color']= [ elem_color[gender] for gender in G.vs['elem']]
	visual_style['vertex_label']= G.vs['name']
	visual_style.update(kwds)
	ig.plot( G, **visual_style )
