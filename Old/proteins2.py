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
