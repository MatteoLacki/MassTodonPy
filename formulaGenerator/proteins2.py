from aminoAcid2 import AminoAcid
import igraph as ig

fasta 	= 'RPKPQQFFGLM'
AAtags 	= ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

# AA = {}
# for aa in AAtags:
# 	AA[aa] = AminoAcid(aa)

x = AminoAcid('A')
print x.superAtoms
x.superAtoms.vs['elementContent']
x.superAtoms.es['Roep']

def elementContent(G):
	atomNo = Counter()	
	for el in G.vs['elem']:
		atomNo[el] += 1
	return atomNo

def fragmentNO(fragment, isProline=False):
	if isProline:
		if 'Calpha' in fragment.vs['name']:
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


def localSuperAtoms(x):
	isProline 	= x.G['name']=='Proline'
	G 			= x.G.copy()
	bonds2break = [ i for i,bt in enumerate(G.es['Roep']) if not bt==None ]
	G.delete_edges(bonds2break)
	if isProline:
		superAtoms = ig.Graph(2)
		superAtoms.add_edge(0,1,Roep='ax')
	else:
		superAtoms = ig.Graph(3)
		superAtoms.add_edge(0,1,Roep='cz')
		superAtoms.add_edge(1,2,Roep='ax')
	for f in x.G.decompose():
		superAtoms.vs[ fragmentNO(f, isProline) ]['elementContent'] = elementContent(f)
	return superAtoms


# W = x.G.copy()

superAtoms.vs['elementContent']

SuperAtoms = localSuperAtoms(x)
SuperAtoms.vs['elementContent']


	# the gist of it
for s in fasta:
	useAppropriateGraph

	if 'by' in bondTypes2severe:


bonds = {}
for i,bondType in enumerate(G.es['Roep']):
	if not bondType == None:
		bonds[bondType] = i

severedBonds = 


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
