import igraph as ig

def style(G, label='name'):
	elem_color = {'C':'grey', 'H':'red', 'N':'blue', 'O':'white', 'S':'pink' }
	bond_thick = {'-':1,'=':8}
	visual_style = {}
	visual_style['vertex_color']= [ elem_color[gender] for gender in G.vs['elem']]
	visual_style['vertex_label']= G.vs[label]
	visual_style['edge_label']= G.es['Roep']
	# visual_style['edge_width'] 	= [ bond_thick[bond] for bond in G.es['bond']]
	return visual_style

def style2(G, label='name'):
	elem_color = {'C':'grey', 'H':'red', 'N':'blue', 'O':'white', 'S':'pink' }
	bond_thick = {'-':1,'=':8}
	visual_style = {}
	visual_style['vertex_color']= [ elem_color[gender] for gender in G.vs['elem']]
	visual_style['vertex_label']= G.vs[label]
	return visual_style

######################################################################

def Backbone( acyclic = True ):
	B = makeAtoms({ 'C':2, 'H':4, 'N':1, 'O':2 })
	B.add_edge( 'N0', 'C0', Roep='cz' )
	B.add_edge( 'C0', 'C1', Roep='ax' )
	B = B + ('N0','H0') + ('N0','H1') + ('C0','H2')+\
	('C1','O0')+('C1','O0')+('C1','O1') + ('O1','H3')	
	B.vs['PDB'] = ['H11', 'H12', 'HA', "H''", 'CA', 'C', "O'", "O''", 'N' ]
	if acyclic:
		Ccarboxyl = 'C1'
		Ocarboxyl1= 'O1_1'
		Ocarboxyl2= 'O1_2'
		Hcarboxyl = 'H1'
	else: 
		Ccarboxyl = 'Ccarboxyl'
		Ocarboxyl1= 'Ocarboxyl1_1'
		Ocarboxyl2= 'Ocarboxyl1_2'
		Hcarboxyl = 'Hcarboxyl1'
	B.vs['IUPAC'] = ['HNalpha1', 'HNalpha2', 'Halpha', Hcarboxyl, 'Calpha', Ccarboxyl, Ocarboxyl1, Ocarboxyl2, 'Nalpha']
	return B

# GG = Backbone( acyclic = False )
# GG.vs['name']
# GG.vs['PDB']
# GG.vs['IUPAC']
# GG.es['Roep']
# ig.plot( GG, **style(GG, label='IUPAC') )
# ig.plot( GG, **style(GG) )
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

# def addBackbone(AA, linkerName):
# 	B = Backbone()
# 	Calpha = B.vs.find('C0').index
# 	linker = AA.vs.find(linkerName).index + len(B.vs)
# 	G = join( B, AA, ['elem'], ['Roep'] )
# 	G = G + ( Calpha, linker )
# 	return G
def addBackbone(A, CalphaNeigbour = 'Cbeta'):
	B = Backbone()
	C = join( A, B, vertex_attributes=['elem', 'IUPAC'], edge_attributes=['Roep'] )
	atomBackbone = C.vs.find(IUPAC='Calpha').index
	atomResidual = C.vs.find(IUPAC=CalphaNeigbour).index
	C.add_edge( source = atomBackbone, target = atomResidual )
	return C

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
	A = makeAtoms({ 'C':1, 'H':3 })+('C0','H0')+('C0','H1')+('C0','H2')
	A.vs['IUPAC'] = ['Hbeta1', 'Hbeta2', 'Hbeta3', 'Cbeta']
	A = addBackbone(A,'Cbeta')
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

R = makeAtoms({ 'C':4, 'H':10, 'N':3 })+\
('C0','H0')+('C0','H1')+('C0','C1')+\
('C1','H2')+('C1','H3')+('C1','C2')+\
('C2','H4')+('C2','H5')+('C2','N0')+\
('N0','H6')+('N0','C3')+\
('C3','N1')+('C3','N1')+('C3','N2')+\
('N1','H7')+\
('N2','H8')+('N2','H9')

R.vs['name']
R.vs['IUPAC'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'Hdelta1', 'Hdelta2', 'HNdelta', 'HNomega1', 'HNomega2_1', 'HNomega2_2', 'Cbeta', 'Cgamma', 'Cdelta', 'guadinino-C', 'Ndelta', 'Nomega1', 'Nomega2']
R = addBackbone(R,'Cbeta')
ig.plot( R, **style2(R,label='IUPAC') )
R['name'] ='Arginine'
return R

# GG = A()
# ig.plot( GG, **style(GG, label='IUPAC') )
# ig.plot( GG, **style(GG) )




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

def I():
	I = makeAtoms({ 'C':4, 'H':9 })+\
	('C0','H0')+('C0','C1')+('C0','C2')+\
	('C1','H1')+('C1','H2')+('C1','H3')+\
	('C2','H4')+('C2','H5')+('C2','C3')+\
	('C3','H6')+('C3','H7')+('C3','H8')
	I = addBackbone(I,'C0')
	I['name'] ='Isoleucine'
	return I

def L():
	L = makeAtoms({ 'C':4, 'H':9 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','H2')+('C1','C2')+('C1','C3')+\
	('C2','H3')+('C2','H4')+('C2','H5')+\
	('C3','H6')+('C3','H7')+('C3','H8')
	L = addBackbone(L,'C0')
	L['name'] ='Leucine'
	return L

def K():
	K = makeAtoms({ 'C':4, 'H':10, 'N':1 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','H2')+('C1','H3')+('C1','C2')+\
	('C2','H4')+('C2','H5')+('C2','C3')+\
	('C3','H6')+('C3','H7')+('C3','N0')+\
	('N0','H8')+('N0','H9')
	K = addBackbone(K,'C0')
	K['name'] ='Lysine'
	return K

def M():
	M = makeAtoms({ 'C':3, 'H':7, 'S':1 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','H2')+('C1','H3')+('C1','S0')+\
	('S0','C2')+\
	('C2','H4')+('C2','H5')+('C2','H6')
	M = addBackbone(M,'C0')
	M['name'] = 'Methionine'
	return M

def F():
	F = makeAtoms({ 'C':7, 'H':7 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','C2')+('C1','C2')+('C1','C6')+\
	('C2','H2')+('C2','C3')+\
	('C3','H3')+('C3','C4')+('C3','C4')+\
	('C4','H4')+('C4','C5')+\
	('C5','H5')+('C5','C6')+('C5','C6')+\
	('C6','H6')
	F = addBackbone(F,'C0')
	F['name'] = 'Phenylalanine'
	return F

def P(): # This is special: links to backbone in two places.
	P = makeAtoms({ 'C':5, 'H':9, 'N':1, 'O':2 })
	P.add_edge( 'C0', 'C1', Roep='ax' )
	P.add_edge( 'N0', 'C1', Roep='cz' )
	P = P + ('N0','H0') + ('C1','H1')+\
	('C0','O0')+('C0','O0')+('C0','O1') + ('O1','H2')+\
	('C1','C2')+\
	('C2','H3')+('C2','H4')+('C2','C3')+\
	('C3','H5')+('C3','H6')+('C3','C4')+\
	('C4','N0')+('C4','H7')+('C4','H8')
	('N')	
	P['name'] = 'Proline'
	return P

def S():
	S = makeAtoms({ 'C':1, 'H':3, 'O':1 })+\
	('C0','H0')+('C0','H1')+('C0','O0')+\
	('O0','H2')
	S = addBackbone(S,'C0')
	S['name'] = 'Serine'
	return S

def T():
	T = makeAtoms({ 'C':2, 'H':5, 'O':1 })+\
	('C0','H0')+('C0','O0')+('C0','C1')+\
	('O0','H1')+\
	('C1','H2')+('C1','H3')+('C1','H4')
	T = addBackbone(T,'C0')
	T['name'] = 'Threonine'
	return T

def W():
	W = makeAtoms({ 'C':9, 'H':8, 'N':1 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','C2')+('C1','C3')+('C1','C3')+\
	('C2','C4')+('C2','C5')+('C2','C5')+\
	('C3','N0')+('C3','H2')+\
	('N0','H3')+('N0','C5')+\
	('C5','C6')+\
	('C6','H4')+('C6','C7')+('C6','C7')+\
	('C4','C8')+('C4','C8')+('C4','H5')+\
	('C7','H6')+('C7','C8')+\
	('C8','H7')
	W = addBackbone(W,'C0')
	W['name'] = 'Tryptophan'
	return W

def Y():
	Y = makeAtoms({ 'C':7, 'H':7, 'O':1 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','C2')+('C1','C2')+('C1','C3')+\
	('C2','H2')+('C2','C4')+\
	('C3','H3')+('C3','C5')+('C3','C5')+\
	('C4','H4')+('C4','C6')+('C4','C6')+\
	('C5','H5')+('C5','C6')+\
	('C6','O0')+\
	('O0','H6')
	Y = addBackbone(Y,'C0')
	Y['name'] = 'Tyrosine'
	return Y

def V():
	V = makeAtoms({ 'C':3, 'H':7 })+\
	('C0','C1')+('C0','C2')+('C0','H0')+\
	('C1','H1')+('C1','H2')+('C1','H3')+\
	('C2','H4')+('C2','H5')+('C2','H6')
	V = addBackbone(V,'C0')
	V['name'] = 'Valine'
	return V


AA = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
# savePath = '/Users/matteo/Dropbox/Science/MassSpectrometry/masstodon/Visual/AminoAcids/'
# savePath = '/Users/matteo/MassSpec/Visual/AminoAcids/'
# for aa in AA:
# 	GG = locals()[aa]()
# 	ig.plot( GG, target = savePath + aa + '.pdf', **style(GG) )