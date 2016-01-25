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

def BondTypes(bonds = ['cz']):
	B2atoms = {	'cz': [ ('Nalpha', 'Calpha'), ('N1', 'C2') ],
				'by': [ ('Nalpha', 'Ccarbo'), ('N1', 'Ccarbo') ],
				'ax': [ ('Calpha','Ccarbo'), ('N1', 'Ccarbo') ]  	}
	B2B = defaultdict(lambda: None)
	try: 
		for b in bonds:
			for s,t in B2atoms[b]:	
				B2B[(s,t)] = B2B[(t,s)] = b
		return B2B
	except KeyError:
		print 'This type of bond is yet not considered: '+str(b)

B2B = BondTypes()

def bonds(fasta, B2B):
	bondTypes= set(B2B[k] for k in B2B)
	notFirst = False
	for left, right, G in aminosIter(fasta):
		for e in G.es:
			yield B2B[ tuple(vName for vName in G.vs[e.tuple]['name']) ]
		if 'by' in bondTypes:
			if notFirst:
				yield 'by'
			else:
				notFirst = True
				yield None
		else:
			yield None

def propGen(propName, fasta):
	if propName == 'name':
		return ( str(aaNo) + ':' + pN for aaNo, (_,_,G) in enumerate( aminosIter(fasta)) for pN in G.vs[propName] )
	else:
		return ( pN for _,_,G in aminosIter(fasta) for pN in G.vs[propName] )

def makeProtein( fasta, attributeNames=('name', 'elem') ):
	N 	= sum( len(G.vs) for l,r,G in aminosIter(fasta) )
	B2B = BondTypes()
	G 	= ig.Graph(N)
	G.add_edges( edges(fasta) )
	for name in attributeNames:
		G.vs[name] = [ a for a in propGen(name, fasta) ]
	G.es['Roep'] = list(bonds(fasta,B2B))
	return G



fasta = substanceP
fasta = 'AA'
G = makeProtein(fasta)		
print G
G.es['Roep']
plott(G, edge_label=G.es['Roep'], bbox=(2000,1000))



for vs in G.vs['name']:

G.vs([0,2,4])['name']

G.vs.select



vertices


G= {'v1': ['v2', 'v3'], 'v2': ['v1'], 'v3': ['v1', 'v4'], 'v4': ['v3']}
mvi= {'v1': 1, 'v2': 2, 'v3': 3, 'v4': 4}
graph= ig.Graph(edges= [(mvi[v], mvi[a]) for v in G.keys() for a in G[v]])
print graph

x.ig.Graph(
	edges 
)

# Try using that.
x = ig.Graph.DictList(
	vertices=[{'name':'A'}, {'name':'B'}], 
	edges 	=[{'s':'A', 't':'B', 'name': 'Hello'}],
	vertex_name_attr = 'name',
	edge_foreign_keys= ('s','t')
)
print x

G.get_eid('9:Ccarbo','10:Nalpha')
print G

for i in G.es[1:10]:
	print i.tuple

for i,s in enumerate(fasta):
	if s=='P':

	else:




# print G
# plott(G, bbox=(1000,1000))

ends, __,starts = G.bfs(0)


