from aminoAcid import AminoAcid
from collections import Counter
import igraph as ig

class SuperAtomGraph(object):
	def __init__(self, fasta='RPKPQQFFGLM'):
		self.fasta 	= fasta
		self.AAtags = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
		self.SuperAtomTypes = {}
		for aa in self.AAtags:
			self.SuperAtomTypes[aa] = self.getSuperAtomType( AminoAcid(aa) )
		self.SuperAtoms = ig.Graph()

	def elementContent(self, G):
		atomNo = Counter()	
		for el in G.vs['elem']:
			atomNo[el] += 1
		return atomNo

	def fragmentNO(self, fragment, isProline=False):
		if isProline:
			if 'C2' in fragment.vs['name']: #Class Amino Acids is not encapsulated
				return 0
			else:
				return 1
		else:
			if 'Calpha' in fragment.vs['name']:
				return 1
			elif 'Nalpha' in fragment.vs['name']:
				return 0
			else:
				return 2

	def getSuperAtomType(self, AA):	
		isProline = AA.G['name']=='Proline'
		G = AA.G.copy()
			# If we add more bonds to break, we will have to modify the aminoAcids.
			# Here we will have to make some intersection with the bond type we wished to break.
		bonds2break = [ i for i,bt in enumerate(G.es['Roep']) if not bt==None ]
		G.delete_edges(bonds2break)
		if isProline:
			superAtoms = ig.Graph(2)
			superAtoms.add_edge(0,1,Roep='ax')
		else:
			superAtoms = ig.Graph(3)
			superAtoms.add_edge(0,1,Roep='cz')
			superAtoms.add_edge(1,2,Roep='ax')
		for f in G.decompose():
			superAtoms.vs[ self.fragmentNO(f, isProline) ]['elementContent'] = self.elementContent(f)
		return superAtoms

	def extendSuperGraph(self, G, newAttr={'Roep':'by'} ):
		N = len(self.SuperAtoms.vs)
		if N==0:
			self.SuperAtoms = G.copy()
		else:
			ends,__,starts = G.bfs(0)
			for v in G.bfsiter(0):
				self.SuperAtoms.add_vertex( **v.attributes() )
				# if len(v.attributes())>0:
				# 	self.SuperAtoms.add_vertex(**v.attributes())
				# else:
				# 	self.SuperAtoms.add_vertex()
			for s,e in zip(starts,ends):
				if e==0:
					self.SuperAtoms.add_edge( source=N-1, target=N, **newAttr )
				else:
					self.SuperAtoms.add_edge( source=N+s, target=N+e, **G.es[ G.get_eid(s,e) ].attributes() )

	def superAtomize(self):
		for s in self.fasta:
			self.extendSuperGraph( self.SuperAtomTypes[s] )

	def plot(self, target=None):
		visual_style = {}
		visual_style['edge_label']	= self.SuperAtoms.es['Roep']
		visual_style['vertex_label']= self.SuperAtoms.vs['elementContent']
		visual_style['bbox'] 		= (1000,2000)
		if target:
			visual_style['target'] = target
		ig.plot( self.SuperAtoms, **visual_style )			

		

		