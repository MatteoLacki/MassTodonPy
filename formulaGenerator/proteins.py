import cPickle as pickle
import igraph as ig

AAS = pickle.load( open('aminoAcidGraphs.txt', 'rb') )

def edges(fasta, AAS):
	prevElemNo = 0
	for s in fasta:
		AA = AAS[s]
		for e in AA.es:
			a, b = e.tuple 
			yield ( a+prevElemNo, b+prevElemNo )
		prevElemNo += len(AA.vs)

def propGen(propName, fasta, AAS):
	if propName == 'name':
		return ( str(aaNo) + '_' + pN for aaNo,s in enumerate(fasta) for pN in AAS[s].vs[propName] )
	else:
		return ( pN for s in fasta for pN in AAS[s].vs[propName] )

def makeProtein(fasta,AAS):
	TotalVno = sum( len(AAS[f].vs) for f in fasta )
	G = ig.Graph(TotalVno)
	G.add_edges( edges(fasta, AAS) )
	for name in ['name', 'elem']:
		G.vs[name] = [ a for a in propGen( name, fasta, AAS ) ]
	return G

G = makeProtein('RPKPQQFFGLM', AAS)		
print G

ig.plot(G)
