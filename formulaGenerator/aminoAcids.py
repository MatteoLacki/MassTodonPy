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
		Hcarboxyl = 'HO1_2'
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
def addBackbone(A, CalphaNeigbour = 'Cbeta', acyclic = True ):
	B = Backbone(acyclic)
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
	R.vs['IUPAC'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'Hdelta1', 'Hdelta2', 'HNdelta', 'HNomega1', 'HNomega2_1', 'HNomega2_2', 'Cbeta', 'Cgamma', 'Cdelta', 'guadinino-C', 'Ndelta', 'Nomega1', 'Nomega2']
	R = addBackbone(R,'Cbeta')
	R['name'] ='Arginine'
	return R

def N():
	N = makeAtoms({ 'C':2, 'H':4, 'N':1, 'O':1 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','O0')+('C1','O0')+('C1','N0')+\
	('N0','H2')+('N0','H3')
	N.vs['IUPAC'] = ['Hbeta1', 'Hbeta2', 'HNgamma1', 'HNgamma2', 'Cbeta', 'Cgamma', 'Ogamma', 'Ngamma']
	N = addBackbone(N,'Cbeta')
	N['name'] ='Asparagine'
	return N

def D():
	D = makeAtoms({ 'C':2, 'H':3, 'O':2 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','O0')+('C1','O0')+('C1','O1')+\
	('O1','H2')
	D.vs['IUPAC'] = ['Hbeta1', 'Hbeta2', 'HOgamma2', 'Cbeta', 'Cgamma', 'Ogamma1', 'Ogamma2']
	D = addBackbone(D,'Cbeta')
	D['name'] ='Aspartic acid'
	return D

def C():
	C = makeAtoms({ 'C':1, 'H':3, 'S':1 })+\
	('C0','H0')+('C0','H1')+('C0','S0')+\
	('S0','H2')
	C.vs['IUPAC'] = ['Hbeta1', 'Hbeta2', 'HSbeta', 'Cbeta', 'Sbeta']
	C = addBackbone(C,'Cbeta')
	C['name'] = 'Cysteine'
	return C

def Q():
	Q = makeAtoms({ 'C':3, 'H':6, 'N':1, 'O':1 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','H2')+('C1','H3')+('C1','C2')+\
	('C2','O0')+('C2','O0')+('C2','N0')+\
	('N0','H4')+('N0','H5')
	Q.vs['IUPAC'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'HNdelta1', 'HNdelta2', 'Cbeta', 'Cgamma', 'Cdelta', 'Odelta', 'Ndelta']
	Q = addBackbone(Q,'Cbeta')
	Q['name'] ='Glutamine'
	return Q

def E():
	E = makeAtoms({ 'C':3, 'H':5, 'O':2 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','H2')+('C1','H3')+('C1','C2')+\
	('C2','O0')+('C2','O0')+('C2','O1')+\
	('O1','H4')
	E.vs['IUPAC'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'HOdelta', 'Cbeta', 'Cgamma', 'Cdelta', 'Odelta1', 'Odelta2']
	E = addBackbone(E,'Cbeta')
	E['name'] ='Glutamic acid'
	return E

def G():
	G = Backbone()
	G.add_vertex('Halpha2',elem='H')
	G.add_edge('Halpha2','C0')
	G.vs['IUPAC'] = ['HNalpha1', 'HNalpha2', 'Halpha1', 'H1', 'Calpha', 'C1', 'O1_1', 'O1_2', 'Nalpha', 'Halpha2']
	return G

def H():
	H = makeAtoms({ 'C':4, 'H':5, 'N':2 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','C2')+('C1','C2')+('C1','N0')+\
	('N0','C3')+('N0','C3')+\
	('C3','H3')+('C3','N1')+('N1','H2')+\
	('N1','C2')+\
	('C2','H4')
	H.vs['IUPAC'] = ['Hbeta1', 'Hbeta2', 'Hpi', 'Htau', 'H4', 'Cbeta', 'C4', 'C5', 'C2', 'Npi', 'Ntau']
	['Hbeta1', 'Hbeta2', 'Htau', 'H3', 'H2', 'Cbeta', 'C4', 'C5', 'C2', 'Npi', 'Ntau']
	H = addBackbone(H,'Cbeta')
	H['name'] ='Histidine'
	return H

def I():
	I = makeAtoms({ 'C':4, 'H':9 })+\
	('C0','H0')+('C0','C1')+('C0','C2')+\
	('C1','H1')+('C1','H2')+('C1','H3')+\
	('C2','H4')+('C2','H5')+('C2','C3')+\
	('C3','H6')+('C3','H7')+('C3','H8')
	I.vs['IUPAC'] = ['Hbeta', 'Hbeta1_1', 'Hbeta1_2', 'Hbeta1_3', 'Hgamma1', 'Hgamma2', 'Hdelta1', 'Hdelta2', 'Hdelta3', 'Cbeta', 'Cbeta1', 'Cgamma', 'Cdelta']
	I = addBackbone(I,'Cbeta')
	I['name'] ='Isoleucine'
	return I

def L():
	L = makeAtoms({ 'C':4, 'H':9 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','H2')+('C1','C2')+('C1','C3')+\
	('C2','H3')+('C2','H4')+('C2','H5')+\
	('C3','H6')+('C3','H7')+('C3','H8')
	L.vs['IUPAC'] = ['Hbeta1', 'Hbeta2', 'Hgamma', 'Hgamma1_1', 'Hgamma1_2', 'Hgamma1_3', 'Hgamma2_1', 'Hgamma2_2', 'Hgamma2_3', 'Cbeta', 'Cgamma', 'Cgamma1', 'Cgamma2']
	L = addBackbone(L,'Cbeta')
	L['name'] ='Leucine'
	return L

def K():
	K = makeAtoms({ 'C':4, 'H':10, 'N':1 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','H2')+('C1','H3')+('C1','C2')+\
	('C2','H4')+('C2','H5')+('C2','C3')+\
	('C3','H6')+('C3','H7')+('C3','N0')+\
	('N0','H8')+('N0','H9')
	K.vs['IUPAC'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'Hdelta1', 'Hdelta2', 'Hepsilon1', 'Hepsilon2', 'HNepsilon1', 'HNepsilon2', 'Cbeta', 'Cgamma', 'Cdelta', 'Cepsilon', 'Nepsilon']
	K = addBackbone(K,'Cbeta')
	K['name'] ='Lysine'
	return K

def M():
	M = makeAtoms({ 'C':3, 'H':7, 'S':1 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','H2')+('C1','H3')+('C1','S0')+\
	('S0','C2')+\
	('C2','H4')+('C2','H5')+('C2','H6')
	M.vs['IUPAC'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'Hepsilon1', 'Hepsilon2', 'Hepsilon3', 'Cbeta', 'Cgamma', 'Cepsilon', 'Sdelta']
	M = addBackbone(M,'Cbeta')
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
	F.vs['IUPAC']=['Hbeta1', 'Hbeta2', 'H6', 'H5', 'H4', 'H3', 'H2', 'Cbeta', 'C1', 'C6', 'C5', 'C4', 'C3', 'C2']
	F = addBackbone(F,'Cbeta',False)
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
	P.vs['IUPAC'] = ['HN1', 'H2', 'HO1_2', 'H3_1', 'H3_2', 'H4_1', 'H4_2', 'H5_1', 'H5_2', 'C1', 'C2', 'C3', 'C4', 'C5', 'O1_1', 'O1_2', 'N1']
	P['name'] = 'Proline'
	return P

def S():
	S = makeAtoms({ 'C':1, 'H':3, 'O':1 })+\
	('C0','H0')+('C0','H1')+('C0','O0')+\
	('O0','H2')
	S.vs['IUPAC']=['Hbeta1', 'Hbeta2', 'HObeta', 'Cbeta', 'Obeta']
	S = addBackbone(S,'Cbeta')
	S['name'] = 'Serine'
	return S

def T():
	T = makeAtoms({ 'C':2, 'H':5, 'O':1 })+\
	('C0','H0')+('C0','O0')+('C0','C1')+\
	('O0','H1')+\
	('C1','H2')+('C1','H3')+('C1','H4')
	T.vs['IUPAC']=['Hbeta', 'HObeta', 'Hgamma1', 'Hgamma2', 'Hgamma3', 'Cbeta', 'Cgamma', 'Obeta']
	T = addBackbone(T,'Cbeta')
	T['name'] = 'Threonine'
	return T

# Alternative tryptophan.
# W = makeAtoms({ 'C':9, 'H':8, 'N':1 })+\
# ('C0','H0')+('C0','H1')+('C0','C1')+\
# ('C1','C2')+('C1','C3')+('C1','C3')+\
# ('C2','C4')+('C2','C5')+('C2','C5')+\
# ('C3','N0')+('C3','H2')+\
# ('N0','H3')+('N0','C5')+\
# ('C5','C6')+\
# ('C6','H4')+('C6','C7')+('C6','C7')+\
# ('C4','C8')+('C4','C8')+('C4','H5')+\
# ('C7','H6')+('C7','C8')+\
# ('C8','H7')
def W():
	W = makeAtoms({ 'C':9, 'H':8, 'N':1 })+\
	('C0','H0')+('C0','H1')+('C0','C1')+\
	('C1','C2')+('C1','C3')+('C1','C3')+\
	('C2','C4')+('C2','C4')+('C2','C5')+\
	('C3','N0')+('C3','H2')+\
	('N0','H3')+('N0','C5')+\
	('C5','C6')+('C5','C6')+\
	('C6','H4')+('C6','C7')+\
	('C4','C8')+('C4','H5')+\
	('C7','H6')+('C7','C8')+('C7','C8')+\
	('C8','H7')
	W.vs['IUPAC']= ['Hbeta1', 'Hbeta2', 'H2', 'HN1', 'H7', 'H4', 'H6', 'H5', 'Cbeta', 'C3', 'C3a', 'C2', 'C4', 'C7a', 'C7', 'C6', 'C5', 'N1']
	W = addBackbone(W,'Cbeta')
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
	Y.vs['IUPAC']=['Hbeta1', 'Hbeta2', 'H6', 'H2', 'H5', 'H3', 'HO4', 'Cbeta', 'C1', 'C6', 'C2', 'C5', 'C3', 'C4', 'O4']
	Y = addBackbone(Y,'Cbeta',False)
	Y['name'] = 'Tyrosine'
	return Y

def V():
	V = makeAtoms({ 'C':3, 'H':7 })+\
	('C0','C1')+('C0','C2')+('C0','H0')+\
	('C1','H1')+('C1','H2')+('C1','H3')+\
	('C2','H4')+('C2','H5')+('C2','H6')
	V.vs['IUPAC']=['Hbeta', 'Hbeta1_1', 'Hbeta1_2', 'Hbeta1_3', 'Hgamma1', 'Hgamma2', 'Hgamma3', 'Cbeta', 'Cbeta1', 'Cgamma']
	V = addBackbone(V,'Cbeta')
	V['name'] = 'Valine'
	return V

AA = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
# savePath = '/Users/matteo/Dropbox/Science/MassSpectrometry/masstodon/Visual/AminoAcids/'
# savePath = '/Users/matteo/MassSpec/Visual/AminoAcids/'
# for aa in AA:
# 	GG = locals()[aa]()
# 	ig.plot( GG, target = savePath + aa + '.pdf', **style(GG) )