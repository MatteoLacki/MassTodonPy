from aminoAcid import AminoAcid
from misc import plott
import igraph as ig

substanceP 	= 'RPKPQQFFGLM'
ubiquitin 	= 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'

def aminosIter(fasta):
	aas 	= ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	AAS 	= {}
	for aa in aas:
		A = AminoAcid(aa)
		AAS[aa] = ( A.getLeft(), A.getRight(), A.getGraph() )
	for i,aa in enumerate(fasta):
		if i==0:
			A = AminoAcid(aa)
			A.add_H()
			yield ( A.getLeft(), A.getRight(), A.getGraph() )
		elif i < len(fasta)-1:
			yield AAS[aa]
		else:
			A = AminoAcid(aa)
			A.add_OH()
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

def propGen(propName, fasta):
	if propName == 'name':
		return ( str(aaNo) + '_' + pN for aaNo, (_,_,G) in enumerate( aminosIter(fasta)) for pN in G.vs[propName] )
	else:
		return ( pN for _,_,G in aminosIter(fasta) for pN in G.vs[propName] )

def makeProtein( fasta, attributeNames=('name', 'elem') ):
	TotalVno = sum( len(G.vs) for l,r,G in aminosIter(fasta) )
	G = ig.Graph(TotalVno)
	G.add_edges( edges(fasta) )
	for name in attributeNames:
		G.vs[name] = [ a for a in propGen(name, fasta) ]
	return G



G = makeProtein()		
plott(G, bbox=(1000,1000))
