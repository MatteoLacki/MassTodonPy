import igraph as ig

def style(G):
	elem_color = {'C':'grey', 'H':'red', 'N':'blue', 'O':'white', 'S':'pink' }
	bond_thick = {'-':1,'=':8}
	visual_style = {}
	visual_style['vertex_color']= [ elem_color[gender] for gender in G.vs['elem']]
	visual_style['vertex_label']= G.vs['elem']
	visual_style['edge_width'] 	= [ bond_thick[bond] for bond in G.es['bond']]
	return visual_style


# def getBackbone():
# 	B = ig.Graph()
# 	B.add_vertex( elem='N', name='N0' )
# 	B.add_vertex( elem='H', name='H0' )
# 	B.add_vertex( elem='H', name='H1' )
# 	B.add_edge( 'N0', 'H0', bond='-' )
# 	B.add_edge( 'N0', 'H1', bond='-' )
# 	B.add_vertex( elem='C', name='C0' ) # C-alpha
# 	B.add_edge( 'C0', 'N0', bond='-' )
# 	B.add_vertex( elem='H', name='H2' )
# 	B.add_edge( 'H2', 'C0', bond='-' )
# 	B.add_vertex( elem='C', name='C1' )
# 	B.add_edge( 'C1', 'C0', bond='-' )
# 	B.add_vertex( elem='O', name='O0' )
# 	B.add_edge( 'O0', 'C1', bond='=' )
# 	B.add_vertex( elem='O', name='O1' )
# 	B.add_edge( 'O1', 'C1', bond='=' )
# 	B.add_vertex( elem='H', name='H3' )
# 	B.add_edge( 'H3', 'O1', bond='-' )

def makeAtoms(atomCnt):

def getBackbone():
	atomCnt = {'C':2, 'H':4, 'N':1, 'O':2}
	vAtt = {'name': [], 'elem': [] }
	for atom in atomCnt:
		vAtt['elem'].extend( [ atom ]*atomCnt[atom] )
		for aC in xrange(atomCnt[atom]):
			vAtt['name'].append( atom+str(aC) )
	B = ig.Graph( sum( atomCnt[a] for a in atomCnt ), vertex_attrs=vAtt )+\
	('N0','H0') + ('N0','H1') + ('C0','N0') + ('C0','H2')+\
	('C0','C1') + ('C1','O0') + ('C1','O1') + ('O1','H3')
	B.es['bond'] = '-'
	B.es[B.get_eid('C1','O0')]['bond']='='
	return B
	
def A():
	A = ig.Graph()
	A.add_vertex( elem='C', name='C0' )
	for i in xrange(3):
		Hname = 'H'+str(i)
		A.add_vertex( elem='H', name=Hname )
		A.add_edge( 'C0', Hname, bond='-' )
	return A

def R():
	R = getBackbone()
	h = 0
	cPrev = 'Calpha'
	for c in xrange(3):
		Cname = 'C'+str(c)
		R.add_vertex( 	elem = 'C', name = Cname )
		R.add_edge( Cname, cPrev )
		for i in xrange(2):
			Hname = 'C'+str(h)
			R.add_vertex( elem = 'H', name = Hname )
			R.add_edge( Cname, Hname )
			h += 1
		cPrev = Cname
	R.add_vertex( elem='N', name='-N(-H)-' )
	R.add_edge( '-N(-H)-', cPrev )

	c += 1
	Cname = 'C'+str(c)
	R.add_vertex( 	elem = 'C', name = Cname )
	R.add_edge( Cname, '-N(-H)-' )
	




# back = getBackbone()
# print back
# back + vert
# ig.plot( back, **style(back) )

# R = getBackbone()
# Hidx = 0
# for i in xrange(3):
# 	Cname = 'C_'+str(i)
# 	R.add_vertex( elem = 'C', 	name = Cname )
# 	R.add_edge( Cname, 'Calpha', bond = '-'  )
# 	for j in xrange(2):
# 		Hname= 'H_'+str(Hidx)
# 		Hidx += 1
# 		R.add_vertex( 	elem='H', 	name=Hname )
# 		R.add_edge( Cname, Hname, bond='-' )


# g = ig.Graph() 
# g.add_vertices(3)
# g.add_edges([(0,2),(1,2)])
# g.vs['name'] = g.vs['type'] = ['a','b','c']

# h = ig.Graph() 
# h.add_vertices(3)
# h.add_edges([(0,1),(1,2)])
# h.vs['name'] = h.vs['type'] = ['a1','b1','c1']
# # We can simply leave out the names. They are not so important. 

# def get_attrs_or_nones(seq, attr_name):
#     try:
#         return seq[attr_name]
#     except KeyError:
#         return [None] * len(seq)

# def better_disjoint_union(g1, g2):
#     g = g1+g2
#     vertex_attributes = set(g.vertex_attributes() + g2.vertex_attributes())
#     edge_attributes = set(g.edge_attributes() + g2.edge_attributes())
#     for attr in vertex_attributes:
#         g.vs[attr] = get_attrs_or_nones(g1.vs, attr) + get_attrs_or_nones(g2.vs, attr)
#     for attr in edge_attributes:
#         g.es[attr] = get_attrs_or_nones(g1.es, attr) + get_attrs_or_nones(g2.es, attr)
#     return g



# print g
# print h
# print better_disjoint_union(g, h)

# testAA = A()
# ig.plot( testAA, **style(testAA) )
# A = ig.Graph()
# A.add_vertex( elem='C' )