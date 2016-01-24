from SuperAtoms import SuperAtomGraph
SAG = SuperAtomGraph(fasta='RPKPQQFFGLM'*100)
SAG.superAtomize()

print SAG.SuperAtoms

SAG.plot()
print SAG.SuperAtoms
SAG.SuperAtoms.vs['elementContent']
SAG.SuperAtoms.es['Roep']


from aminoAcid import AminoAcid
import igraph as ig

x = AminoAcid('A')
G = x.G.copy()

# Write it so that it follows the order.
N = ig.Graph(len(G.vs)*1000)
for i in xrange(1000):


len(G.vs)
for i in xrange(1000):
	G += G




EC = SAG.SuperAtoms.vs['elementContent']



SAG.SuperAtomTypes['P'].vs['elementContent']


def elementContent( G):
	atomNo = Counter()	
	for el in G.vs['elem']:
		atomNo[el] += 1
	return atomNo

def fragmentNO(fragment, isProline=False):
	if isProline:
		if 'C2' in fragment.vs['name']:
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

AA = AminoAcid('P')
def getSuperAtomType(self, AA):	
	isProline = AA.G['name']=='Proline'
	G = AA.G.copy()
		# If we add more bonds to break, we will have to modify the aminoAcids.
		# Here we will have to make some intersection with the bond type we wished to break.
G.es['Roep']
	bonds2break = [ i for i,bt in enumerate(G.es['Roep']) if not bt==None ]
	G.delete_edges(bonds2break)
print G
ig.plot(G)
	if isProline:
		superAtoms = ig.Graph(2)
		superAtoms.add_edge(0,1,Roep='ax')
print superAtoms
	else:
		superAtoms = ig.Graph(3)
		superAtoms.add_edge(0,1,Roep='cz')
		superAtoms.add_edge(1,2,Roep='ax')
	for f in G.decompose():
		superAtoms.vs[ fragmentNO(f, isProline) ]['elementContent'] = elementContent(f)
	return superAtoms
superAtoms.vs['elementContent']
y = G.decompose()
fragment = y[0]

from aminoAcid import AminoAcid






def extendG(G, E, newAttr={'name':'link','col':'G'} ):
	N = len(G.vs)
	if N==0:
		return E.copy()
	else:
		ends, __,starts  = E.bfs(0)
		for v in E.bfsiter(0):
			if len(v.attributes())>0:
				G.add_vertex(**v.attributes())
			else:
				G.add_vertex()
		for s,e in zip(starts,ends):
			if e==0:
				G.add_edge( source=N-1, target=N, **newAttr )
			else:
				G.add_edge( source=N+s, target=N+e, **E.es[ E.get_eid(s,e) ].attributes() )
		return G

x = AminoAcid('A')
y = AminoAcid('R')
z = AminoAcid('P')

superAtoms = ig.Graph()
print superAtoms
superAtoms = extendG(superAtoms, x.superAtoms, {'Roep':'by'})
print superAtoms
superAtoms = extendG(superAtoms, y.superAtoms, {'Roep':'by'})
print superAtoms
superAtoms = extendG(superAtoms, z.superAtoms, {'Roep':'by'})
print superAtoms
superAtoms.vs['elementContent']


	# the gist of it
for s in fasta:
	useAppropriateGraph

	if 'by' in bondTypes2severe:


class superAtoms(object):
	def __init__(self, bondTypes2severe): 
		self.G = ig.Graph()
		self.lastSuperAtom = None
		self.bondTypes2severe = set(bondTypes2severe)
	def addAA(self, AA):



# for aa in aas:
# 	AminoAcid(aa).plot(target='/Volumes/doom/Users/matteo/Dropbox/Science/MassSpectrometry/MassTodon/Visual/AAs/'+aa+'.pdf')

# try:
# 	for aa in aas:
# 		X = AminoAcid(aa)
# 		X.add_OH()
# 		X.add_H()
# 		X.plot(target='/Volumes/doom/Users/matteo/Dropbox/Science/MassSpectrometry/MassTodon/Visual/AAs/'+aa+'_with_OH_and_H.pdf')
# except ValueError:
# 	print aa
