import igraph as ig

def style(G):
	elem_color = {'C':'grey', 'H':'red', 'N':'blue', 'O':'white', 'S':'pink' }
	bond_thick = {'-':1,'=':8}
	visual_style = {}
	visual_style['vertex_color']= [ elem_color[gender] for gender in G.vs['elem']]
	visual_style['vertex_label']= G.vs['elem']
	# visual_style['edge_width'] 	= [ bond_thick[bond] for bond in G.es['bond']]
	return visual_style

######################################################################

def Backbone():
	B = makeAtoms({ 'C':2, 'H':4, 'N':1, 'O':2 })
	B.add_edge( 'N0', 'C0', Roep='cz' )
	B.add_edge( 'C0', 'C1', Roep='ax' )
	B = B + ('N0','H0') + ('N0','H1') + ('C0','H2')+\
	('C1','O0')+('C1','O0')+('C1','O1') + ('O1','H3')	
	return B

######################################################################

def getAttr(seq, attr_name):
    try:
        return seq[attr_name]
    except KeyError:
        return [None] * len(seq)

def join(g1, g2, vertex_attributes, edge_attributes):
    g = g1+g2
    for attr in vertex_attributes:
        g.vs[attr] = getAttr(g1.vs, attr) + getAttr(g2.vs, attr)
    for attr in edge_attributes:
        g.es[attr] = getAttr(g1.es, attr) + getAttr(g2.es, attr)
    return g

def addBackbone(AA, linkerName):
	B = Backbone()
	Calpha = B.vs.find('C0').index
	linker = AA.vs.find(linkerName).index + len(B.vs)
	G = join( B, AA, ['elem'], ['Roep'] )
	G = G + ( Calpha, linker )
	return G

def makeAtoms(atomCnt):
	vAtt = {'name': [], 'elem': [] }
	for atom in atomCnt:
		vAtt['elem'].extend( [ atom ]*atomCnt[atom] )
		for aC in xrange(atomCnt[atom]):
			vAtt['name'].append( atom+str(aC) )
	G = ig.Graph( sum( atomCnt[a] for a in atomCnt ), vertex_attrs=vAtt )
	return G
######################################################################

def A():
	A = makeAtoms({ 'C':1, 'H':3 })+\
	('C0','H0')+('C0','H1')+('C0','H2')
	A = addBackbone(A,'C0')
	A['name'] = 'Alanine'
	return A

def R():
	R = makeAtoms({ 'C':4, 'H':10, 'N':3 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','H2')+('C1','H3')+('C1','C2')+\
	('C2','H4')+('C2','H5')+('C2','N0')+\
	('N0','H6')+('N0','C3')+\
	('C3','N1')+('C3','N1')+('C3','N2')+\
	('N1','H7')+\
	('N2','H8')+('N2','H9')
	R = addBackbone(R,'C0')
	R['name'] ='Arginine'
	return R

def N():
	N = makeAtoms({ 'C':2, 'H':4, 'N':1, 'O':1 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','O0')+('C1','O0')+('C1','N0')+\
	('N0','H2')+('N0','H3')
	N = addBackbone(N,'C0')
	N['name'] ='Asparagine'
	return N
def D():
	D = makeAtoms({ 'C':2, 'H':3, 'O':2 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','O0')+('C1','O0')+('C1','O1')+\
	('O1','H2')
	D = addBackbone(D,'C0')
	D['name'] ='Aspartic acid'
	return D

def C():
	C = makeAtoms({ 'C':1, 'H':3, 'S':1 })+\
	('C0','H0')+('C0','H1')+('C0','S0')+\
	('S0','H2')
	C = addBackbone(C,'C0')
	C['name'] = 'Cysteine'
	return C

def Q():
	Q = makeAtoms({ 'C':3, 'H':6, 'N':1, 'O':1 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','H2')+('C1','H3')+('C1','C2')+\
	('C2','O0')+('C2','O0')+('C2','N0')+\
	('N0','H4')+('N0','H5')
	Q = addBackbone(Q,'C0')
	Q['name'] ='Glutamine'
	return Q

def E():
	E = makeAtoms({ 'C':3, 'H':5, 'O':2 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','H2')+('C1','H3')+('C1','C2')+\
	('C2','O0')+('C2','O0')+('C2','O1')+\
	('O1','H4')
	E = addBackbone(E,'C0')
	E['name'] ='Glutamic acid'
	return E

def G():
	G = makeAtoms({ 'H':1 })
	G = addBackbone(G,'H0')
	G['name'] = 'Glycine'
	return G


def H():
	H = makeAtoms({ 'C':4, 'H':5, 'N':2 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','C2')+('C1','C2')+('C1','N0')+\
	('N0','H2')+('N0','C3')+\
	('C3','H3')+('C3','N1')+('C3','N1')+\
	('N1','C2')+\
	('C2','H4')
	H = addBackbone(H,'C0')
	H['name'] ='Histidine'
	return H

# Add Isoleucine

GG = H()
# { 'C':2, 'H':2, 'N':1, 'O':1 } # The backbone atom count
# GG = H()
print GG 
ig.plot( GG, **style(GG) )
# locals()["myfunction"]()