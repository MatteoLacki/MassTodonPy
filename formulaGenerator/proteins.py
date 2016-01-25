from aminoAcid import AminoAcid
from misc import plott
import igraph as ig
from collections import defaultdict

substanceP 	= 'RPKPQQFFGLM'
ubiquitin 	= 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'

def aminosIter(fasta):
	aas 	= ('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')
	AAS 	= {}
	for aa in aas:
		AAS[aa] = AminoAcid(aa)
	for i,aa in enumerate(fasta):
		if i in (0, len(fasta)-1):
			A = AminoAcid(aa)
			if i == 0:
				A.add_H()
			if i == len(fasta)-1:
				A.add_OH()
		else:
			A = AAS[aa]
		yield ( A.getLeft(), A.getRight(), A.getGraph() )

def edges(fasta):
	prevElemNo = 0
	for left, right, G in aminosIter(fasta):
		for e in G.es:
			a, b = e.tuple 
			yield ( a+prevElemNo, b+prevElemNo ) # bonds within new AA
		if prevElemNo:
			yield ( prevElemNo+left, prevRight ) # bond between last AA and new AA
		prevRight 	= right + prevElemNo
		prevElemNo += len(G.vs)

def BondTypes(bonds = ['cz']):
	B2atoms = {	'cz': [ ('Nalpha', 'Calpha'), ('N1', 'C2') ],
				'by': [ ('Nalpha', 'Ccarbo'), ('N1', 'Ccarbo') ],
				'ax': [ ('Calpha','Ccarbo'), ('C2', 'Ccarbo') ]  	}
	B2B = defaultdict(lambda: None)
	try: 
		for b in bonds:
			for s,t in B2atoms[b]:	
				B2B[(s,t)] = B2B[(t,s)] = b
		return B2B
	except KeyError:
		print 'This type of bond is yet not considered: '+str(b)

def GetBonds(fasta, B2B):
	bondTypes= set(B2B[k] for k in B2B)
	notFirst = False
	for left, right, G in aminosIter(fasta):
		for e in G.es:
			yield B2B[ tuple(vName for vName in G.vs[e.tuple]['name']) ]
		if notFirst:
			if 'by' in bondTypes:
				yield 'by'
			else:
				yield None			
		else:	
			notFirst = True

def propGen(propName, fasta):
	if propName == 'name':
		return ( str(aaNo) + ':' + pN for aaNo, (_,_,G) in enumerate( aminosIter(fasta)) for pN in G.vs[propName] )
	else:
		return ( pN for _,_,G in aminosIter(fasta) for pN in G.vs[propName] )

def makeProtein( fasta, attributeNames=('name', 'elem'), bonds=['cz','ax','by'] ):
	N 	= sum( len(G.vs) for l,r,G in aminosIter(fasta) )
	B2B = BondTypes(bonds)
	G 	= ig.Graph(N)
	G.add_edges( edges(fasta) )
	for name in attributeNames:
		G.vs[name] = [ a for a in propGen(name, fasta) ]
	G.es['Roep'] = list(GetBonds(fasta,B2B))
	return G


fasta 	= substanceP
G 		= makeProtein(fasta, bonds=['cz','by'])		
plott(G, edge_label=G.es['Roep'], bbox=(2000,1000))


# ends, __,starts = G.bfs(0)


