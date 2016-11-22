import igraph as ig

class MissingAminoAcid(Exception): pass
class WrongArgument(Exception): pass
class OH_already_deleted(Exception): pass

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

class AminoAcid(object):
	def __init__(self, aa):
		try:
			self.G = getattr(self, aa)()
		except AttributeError:
			print aa + ' is not among acceptable amino acids:'
			print ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
			raise MissingAminoAcid

		if aa	== 'F' or aa == 'Y': # the only amino acids with ambigous IUPAC notation
			self.OH = ( 'Ocarboxyl1_2', 'Hcarboxyl1')
			self.C1 = 'Ccarboxyl'
		else:
			self.OH = ( 'O1_2', 'HO1_2')
			self.C1 = 'C1'

			# Now we rename all the things.
		self.G.vs['IUPAC'] = self.G.vs['name']
		self.updateVertexNames(0)
		self.AAnoNext = 1

	def __str__(self):
		return self.G.__str__()

	def __repr__(self):
		return self.G.__repr__()

	def addBackbone(self, A, acyclic=True):
		B = self.Backbone(acyclic)
		C = join( A, B, vertex_attributes=['elem','name'], edge_attributes=['Roep'] )
		C += ('Calpha', 'Cbeta')
		return C

	def makeAtoms(self, atomCnt):
		vAtt = {'name': [], 'elem': [] }
		for atom in atomCnt:
			vAtt['elem'].extend( [ atom ]*atomCnt[atom] )
			for aC in xrange(atomCnt[atom]):
				vAtt['name'].append( atom+str(aC) )
		G = ig.Graph( sum( atomCnt[a] for a in atomCnt ), vertex_attrs=vAtt )
		return G

	def Backbone(self, acyclic=True):
		B = self.makeAtoms({ 'C':2, 'H':4, 'N':1, 'O':2 })
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
		B.vs['name'] = ['HNalpha1', 'HNalpha2', 'Halpha', Hcarboxyl, 'Calpha', Ccarboxyl, Ocarboxyl1, Ocarboxyl2, 'Nalpha']
		return B

	def A(self):
		A = self.makeAtoms({'C':1, 'H':3})+('C0','H0')+('C0','H1')+('C0','H2')
		A.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hbeta3', 'Cbeta']
		A = self.addBackbone(A)
		A['name'] = 'Alanine'
		return A

	def R(self):
		R = self.makeAtoms({ 'C':4, 'H':10, 'N':3 })+\
		('C0','H0')+('C0','H1')+('C0','C1')+\
		('C1','H2')+('C1','H3')+('C1','C2')+\
		('C2','H4')+('C2','H5')+('C2','N0')+\
		('N0','H6')+('N0','C3')+\
		('C3','N1')+('C3','N1')+('C3','N2')+\
		('N1','H7')+\
		('N2','H8')+('N2','H9')
		R.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'Hdelta1', 'Hdelta2', 'HNdelta', 'HNomega1', 'HNomega2_1', 'HNomega2_2', 'Cbeta', 'Cgamma', 'Cdelta', 'guadinino-C', 'Ndelta', 'Nomega1', 'Nomega2']
		R = self.addBackbone(R)
		R['name'] ='Arginine'
		return R

	def N(self):
		N = self.makeAtoms({ 'C':2, 'H':4, 'N':1, 'O':1 })+\
		('C0','H0')+('C0','H1')+('C0','C1')+\
		('C1','O0')+('C1','O0')+('C1','N0')+\
		('N0','H2')+('N0','H3')
		N.vs['name'] = ['Hbeta1', 'Hbeta2', 'HNgamma1', 'HNgamma2', 'Cbeta', 'Cgamma', 'Ogamma', 'Ngamma']
		N = self.addBackbone(N)
		N['name'] ='Asparagine'
		return N

	def D(self):
		D = self.makeAtoms({ 'C':2, 'H':3, 'O':2 })+\
		('C0','H0')+('C0','H1')+('C0','C1')+\
		('C1','O0')+('C1','O0')+('C1','O1')+\
		('O1','H2')
		D.vs['name'] = ['Hbeta1', 'Hbeta2', 'HOgamma2', 'Cbeta', 'Cgamma', 'Ogamma1', 'Ogamma2']
		D = self.addBackbone(D)
		D['name'] ='Aspartic acid'
		return D

	def C(self):
		C = self.makeAtoms({ 'C':1, 'H':3, 'S':1 })+\
		('C0','H0')+('C0','H1')+('C0','S0')+\
		('S0','H2')
		C.vs['name'] = ['Hbeta1', 'Hbeta2', 'HSbeta', 'Cbeta', 'Sbeta']
		C = self.addBackbone(C)
		C['name'] = 'Cysteine'
		return C

	def Q(self):
		Q = self.makeAtoms({ 'C':3, 'H':6, 'N':1, 'O':1 })+\
		('C0','H0')+('C0','H1')+('C0','C1')+\
		('C1','H2')+('C1','H3')+('C1','C2')+\
		('C2','O0')+('C2','O0')+('C2','N0')+\
		('N0','H4')+('N0','H5')
		Q.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'HNdelta1', 'HNdelta2', 'Cbeta', 'Cgamma', 'Cdelta', 'Odelta', 'Ndelta']
		Q = self.addBackbone(Q)
		Q['name'] ='Glutamine'
		return Q

	def E(self):
		E = self.makeAtoms({ 'C':3, 'H':5, 'O':2 })+\
		('C0','H0')+('C0','H1')+('C0','C1')+\
		('C1','H2')+('C1','H3')+('C1','C2')+\
		('C2','O0')+('C2','O0')+('C2','O1')+\
		('O1','H4')
		E.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'HOdelta', 'Cbeta', 'Cgamma', 'Cdelta', 'Odelta1', 'Odelta2']
		E = self.addBackbone(E)
		E['name'] ='Glutamic acid'
		return E

	def G(self):
		G = self.Backbone()
		G.add_vertex('Halpha2',elem='H')
		G.add_edge('Halpha2','Calpha')
		G.vs['name'] = ['HNalpha1', 'HNalpha2', 'Halpha1', 'H1', 'Calpha', 'C1', 'O1_1', 'O1_2', 'Nalpha', 'Halpha2']
		G['name'] ='Glycine'
		return G

	def H(self):
		H = self.makeAtoms({ 'C':4, 'H':5, 'N':2 })+\
		('C0','H0')+('C0','H1')+('C0','C1')+\
		('C1','C2')+('C1','C2')+('C1','N0')+\
		('N0','C3')+('N0','C3')+\
		('C3','H3')+('C3','N1')+('N1','H2')+\
		('N1','C2')+\
		('C2','H4')
		H.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hpi', 'Htau', 'H4', 'Cbeta', 'C4', 'C5', 'C2', 'Npi', 'Ntau']
		['Hbeta1', 'Hbeta2', 'Htau', 'H3', 'H2', 'Cbeta', 'C4', 'C5', 'C2', 'Npi', 'Ntau']
		H = self.addBackbone(H)
		H['name'] ='Histidine'
		return H

	def I(self):
		I = self.makeAtoms({ 'C':4, 'H':9 })+\
		('C0','H0')+('C0','C1')+('C0','C2')+\
		('C1','H1')+('C1','H2')+('C1','H3')+\
		('C2','H4')+('C2','H5')+('C2','C3')+\
		('C3','H6')+('C3','H7')+('C3','H8')
		I.vs['name'] = ['Hbeta', 'Hbeta1_1', 'Hbeta1_2', 'Hbeta1_3', 'Hgamma1', 'Hgamma2', 'Hdelta1', 'Hdelta2', 'Hdelta3', 'Cbeta', 'Cbeta1', 'Cgamma', 'Cdelta']
		I = self.addBackbone(I)
		I['name'] ='Isoleucine'
		return I

	def L(self):
		L = self.makeAtoms({ 'C':4, 'H':9 })+\
		('C0','H0')+('C0','H1')+('C0','C1')+\
		('C1','H2')+('C1','C2')+('C1','C3')+\
		('C2','H3')+('C2','H4')+('C2','H5')+\
		('C3','H6')+('C3','H7')+('C3','H8')
		L.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hgamma', 'Hgamma1_1', 'Hgamma1_2', 'Hgamma1_3', 'Hgamma2_1', 'Hgamma2_2', 'Hgamma2_3', 'Cbeta', 'Cgamma', 'Cgamma1', 'Cgamma2']
		L = self.addBackbone(L)
		L['name'] ='Leucine'
		return L

	def K(self):
		K = self.makeAtoms({ 'C':4, 'H':10, 'N':1 })+\
		('C0','H0')+('C0','H1')+('C0','C1')+\
		('C1','H2')+('C1','H3')+('C1','C2')+\
		('C2','H4')+('C2','H5')+('C2','C3')+\
		('C3','H6')+('C3','H7')+('C3','N0')+\
		('N0','H8')+('N0','H9')
		K.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'Hdelta1', 'Hdelta2', 'Hepsilon1', 'Hepsilon2', 'HNepsilon1', 'HNepsilon2', 'Cbeta', 'Cgamma', 'Cdelta', 'Cepsilon', 'Nepsilon']
		K = self.addBackbone(K)
		K['name'] ='Lysine'
		return K

	def M(self):
		M = self.makeAtoms({ 'C':3, 'H':7, 'S':1 })+\
		('C0','H0')+('C0','H1')+('C0','C1')+\
		('C1','H2')+('C1','H3')+('C1','S0')+\
		('S0','C2')+\
		('C2','H4')+('C2','H5')+('C2','H6')
		M.vs['name'] = ['Hbeta1', 'Hbeta2', 'Hgamma1', 'Hgamma2', 'Hepsilon1', 'Hepsilon2', 'Hepsilon3', 'Cbeta', 'Cgamma', 'Cepsilon', 'Sdelta']
		M = self.addBackbone(M)
		M['name'] = 'Methionine'
		return M

	def F(self):
		F = self.makeAtoms({ 'C':7, 'H':7 })+\
		('C0','H0')+('C0','H1')+('C0','C1')+\
		('C1','C2')+('C1','C2')+('C1','C6')+\
		('C2','H2')+('C2','C3')+\
		('C3','H3')+('C3','C4')+('C3','C4')+\
		('C4','H4')+('C4','C5')+\
		('C5','H5')+('C5','C6')+('C5','C6')+\
		('C6','H6')
		F.vs['name']=['Hbeta1', 'Hbeta2', 'H6', 'H5', 'H4', 'H3', 'H2', 'Cbeta', 'C1', 'C6', 'C5', 'C4', 'C3', 'C2']
		F = self.addBackbone(F,False)
		F['name'] = 'Phenylalanine'
		return F

	def P(self): # This is special: links to backbone in two places.
		P = self.makeAtoms({ 'C':5, 'H':9, 'N':1, 'O':2 })
		P.add_edge( 'C0', 'C1', Roep='ax' )
		P.add_edge( 'N0', 'C1', Roep='cz' )
		P = P + ('N0','H0') + ('C1','H1')+\
		('C0','O0')+('C0','O0')+('C0','O1') + ('O1','H2')+\
		('C1','C2')+\
		('C2','H3')+('C2','H4')+('C2','C3')+\
		('C3','H5')+('C3','H6')+('C3','C4')+\
		('C4','N0')+('C4','H7')+('C4','H8')
		P.vs['name'] = ['HN1', 'H2', 'HO1_2', 'H3_1', 'H3_2', 'H4_1', 'H4_2', 'H5_1', 'H5_2', 'C1', 'C2', 'C3', 'C4', 'C5', 'O1_1', 'O1_2', 'N1']
		P['name'] = 'Proline'
		return P

	def S(self):
		S = self.makeAtoms({ 'C':1, 'H':3, 'O':1 })+\
		('C0','H0')+('C0','H1')+('C0','O0')+\
		('O0','H2')
		S.vs['name']=['Hbeta1', 'Hbeta2', 'HObeta', 'Cbeta', 'Obeta']
		S = self.addBackbone(S)
		S['name'] = 'Serine'
		return S

	def T(self):
		T = self.makeAtoms({ 'C':2, 'H':5, 'O':1 })+\
		('C0','H0')+('C0','O0')+('C0','C1')+\
		('O0','H1')+\
		('C1','H2')+('C1','H3')+('C1','H4')
		T.vs['name']=['Hbeta', 'HObeta', 'Hgamma1', 'Hgamma2', 'Hgamma3', 'Cbeta', 'Cgamma', 'Obeta']
		T = self.addBackbone(T)
		T['name'] = 'Threonine'
		return T

	# Alternative tryptophan.
	# W = self.makeAtoms({ 'C':9, 'H':8, 'N':1 })+\
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
	def W(self):
		W = self.makeAtoms({ 'C':9, 'H':8, 'N':1 })+\
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
		W.vs['name']= ['Hbeta1', 'Hbeta2', 'H2', 'HN1', 'H7', 'H4', 'H6', 'H5', 'Cbeta', 'C3', 'C3a', 'C2', 'C4', 'C7a', 'C7', 'C6', 'C5', 'N1']
		W = self.addBackbone(W)
		W['name'] = 'Tryptophan'
		return W

	def Y(self):
		Y = self.makeAtoms({ 'C':7, 'H':7, 'O':1 })+\
		('C0','H0')+('C0','H1')+('C0','C1')+\
		('C1','C2')+('C1','C2')+('C1','C3')+\
		('C2','H2')+('C2','C4')+\
		('C3','H3')+('C3','C5')+('C3','C5')+\
		('C4','H4')+('C4','C6')+('C4','C6')+\
		('C5','H5')+('C5','C6')+\
		('C6','O0')+\
		('O0','H6')
		Y.vs['name']=['Hbeta1', 'Hbeta2', 'H6', 'H2', 'H5', 'H3', 'HO4', 'Cbeta', 'C1', 'C6', 'C2', 'C5', 'C3', 'C4', 'O4']
		Y = self.addBackbone(Y,False)
		Y['name'] = 'Tyrosine'
		return Y

	def V(self):
		V = self.makeAtoms({ 'C':3, 'H':7 })+\
		('C0','C1')+('C0','C2')+('C0','H0')+\
		('C1','H1')+('C1','H2')+('C1','H3')+\
		('C2','H4')+('C2','H5')+('C2','H6')
		V.vs['name']=['Hbeta', 'Hbeta1_1', 'Hbeta1_2', 'Hbeta1_3', 'Hgamma1', 'Hgamma2', 'Hgamma3', 'Cbeta', 'Cbeta1', 'Cgamma']
		V = self.addBackbone(V)
		V['name'] = 'Valine'
		return V

	def plot(self):
		elem_color = {'C':'grey', 'H':'red', 'N':'blue', 'O':'white', 'S':'pink' }
		visual_style = {}
		visual_style['vertex_color']= [ elem_color[gender] for gender in self.G.vs['elem']]
		visual_style['vertex_label']= self.G.vs['name']
		visual_style['edge_label']	= self.G.es['Roep']
		ig.plot( self.G, **visual_style )

	def delOH(self):
		AAno = str(self.AAnoNext-1)
		self.G.delete_vertices([ AAno+'_'+name for name in self.OH ])

	def delH(self):
		if self.AAnoNext > 1:
			raise OH_already_deleted

		# AAno = str(self.AAnoNext-1)
		if self.G['name'] == 'Proline':
			self.G.delete_vertices( self.G.vs.find( '0_HN1' ).index )
		else:
			self.G.delete_vertices( self.G.vs.find( '0_HNalpha1' ).index )
			self.G.vs.find( '0_HNalpha2' )['IUPAC'] = 'HNalpha'

	def updateVertexNames(self, AAno):
		self.AAno = AAno
		self.G.vs['name'] = [ str(AAno)+'_'+name for name in self.G.vs['IUPAC'] ]

	def __iadd__(self, RightAA):
		if isinstance(RightAA, AminoAcid):
			RightAA.delH()
			self.delOH()
			RightAA.updateVertexNames(self.AAnoNext)

			self.G = join( self.G, RightAA.G,
				vertex_attributes	= ['elem','name','IUPAC'],
				edge_attributes		= ['Roep'] )

			AAno = str(self.AAnoNext-1)
			AAnext = str(self.AAnoNext)

			if RightAA.G['name'] == 'Proline':
				self.G.add_edge( AAno+'_'+self.C1, AAnext+'_N1', Roep='by' )
			else:
				self.G.add_edge( AAno+'_'+self.C1, AAnext+'_Nalpha', Roep='by' )

				# Updating the IUPAC-like names of OH.
			self.OH = RightAA.OH
			self.AAnoNext += 1
			return self
		else:
			print 'Ye cannot have +(Object,'+str(type(RightAA))+')'
			raise WrongArgument
