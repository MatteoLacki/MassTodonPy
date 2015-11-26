import intervaltree as it
from exampleData import theory, tol, minExpSupp, empiria, savePath
from collections import Counter, defaultdict
import igraph as ig

BFG = ig.Graph() 					# the BIG FUCKING GRAPH
tolIntervals = it.IntervalTree()

for famNo, spectrum in enumerate(theory):
	famName = 'F_' + str(famNo)
	BFG.add_vertex( type='F', name=famName, potentialSupport=0.0 )
	for peakNo, info in enumerate(spectrum):
		mass, intensity = info
		peakName = 'P_' + str(peakNo) + '_' + str(famNo)
		BFG.add_vertex( type='P',name=peakName, mass=mass, intensity=intensity )
		BFG.add_edge( famName, peakName )
		tolIntervals[ mass-tol : mass+tol ] = peakName

for expPeakNo, info in enumerate(empiria):
	mass, intensity = info
	expPeakName = 'E_' + str(expPeakNo)
	BFG.add_vertex( type='E', name=expPeakName, mass=mass, intensity=intensity )
	for interval in tolIntervals[mass]:
		peakName = interval.data
		BFG.add_edge( expPeakName, peakName )

def getAtheoreticalPeaks(G):
	'Masses and intensities of peaks that cannot be attributed to any theoretically observable chemical compound.'
	return [ (experim['mass'], experim['intensity'] ) for experim in BFG.vs(type='E') ]

def getUnsupportedFamilies(G, minimalProbability):
	'Info on the chemical compounds with potentially less support in data than the minimalProbability value.'
	print 'IMPLEMENT IT'
	pass

for family in BFG.vs(type='F'): 							# Establishing experimental support for individual compounds.	
	for peakNo in BFG.neighbors(family):
		if BFG.neighborhood_size( peakNo, order=1 ) > 2: 	# > 2 for the PEAK and FAMILY peak are within NEIGBORHOOD( ORDER = 1 )
			family['potentialSupport'] += BFG.vs[peakNo]['intensity'] 	
	if family['potentialSupport'] < minExpSupp:
		for peakNo in BFG.neighbors(family):
	 		BFG.delete_edges( 
	 			BFG.get_eid(peakNo, nodeNo) for nodeNo in BFG.neighbors( peakNo ) 
	 			if nodeNo != family.index ) 			# Severe edges between peaks and experimental peaks, not between peaks and family node

subproblems = [ G for G in BFG.decompose() 
				if {'F','E'} <= set(G.vs['type']) ]		# Only good subproblems have both a family node and an empirical node 

G = subproblems[1]
G.es['varNo'] = None

G.add_vertex( type='R',name='Root') 	# Add the root: its edges store the mixture variables (alphas) in the optimisation problem
RootNo = G.vs.find(name='Root').index
for varNo, f in enumerate(G.vs( type_eq='F' )):
	G.add_edge( 'Root', f['name'], varNo=varNo, varType='alpha' ) 	# Number the alphas
varNoMax = varNo+1

AeqSparse = []	# Prepare a scarse matrix representation of the equality constraints in the optimisation problem
for i,P in enumerate( G.vs( type = 'P' ) ):
	for nodeNo in G.neighbors(P):
		node = G.vs[nodeNo]
		isFamilyNode = node['type'] == 'F'
		if isFamilyNode:
			edge = G.es[ G.get_eid(nodeNo, RootNo) ]
		else:	
			edge = G.es[ G.get_eid( P.index, nodeNo) ]
		
		if edge['varNo'] == None: 			
			edge['varNo'] 	= varNoMax 	# Assign number to the flow variable in the optimisation problem
			edge['varType'] = 'flow' 	# Tag as flow
			varNoMax += 1
		
		if isFamilyNode:
			AeqSparse.append(( i, edge['varNo'], -float(P['intensity']) ))
		else:
			AeqSparse.append(( i, edge['varNo'], 1.0 ))

print AeqSparse	# TO DO check, if this is consistent with the tree numeration scheme. Plot tree with numbers on edges 

AineqSparse = [] 	# Prepare a scarse matrix representation of the inequality constraints in the optimisation problem
for i, E in enumerate( G.vs( type = 'E' ) ):
	for e in G.es[ G.incident(E) ]:
		AineqSparse.append(( i, e['varNo'], 1.0 ))
print AineqSparse	# TO DO check, if this is consistent with the tree numeration scheme. Plot tree with numbers on edges 	

cSparse = [ (e['varNo'], 1.) for e in G.es( varType = 'flow' ) ]
print cSparse 	# TO DO check, if this is consistent with the tree numeration scheme
	
# TO DO make dense matrices and plug them into solver



# layout = BFG.layout("kk") 				# Makes a full graph
# G.vs['label'] = G.vs['type']
# color_dict = {'F': 'green', 'P': 'blue', 'E': 'red', 'R': 'yellow' }
# G.vs['color'] = [ color_dict[ge] for ge in G.vs['type']]
# ig.rglplot(G, layout = layout, target=savePath+'/artificialNode.pdf')#,  bbox = (300, 300), margin = 20)


# for i,p in enumerate( G.vs( type_eq='P' ) ):
# 	print p
	# for e in G.es(G.incident(p)):
	# 	print e.tuple
	# print
	# for e in G.incident(p):
		
	# for n in G.neighbors(p):
	# 	print (p.index, n)
	# 	print G.es[n].tuple
	# 	print
		# for h in G.es(p.index, n)
		# 	print h

# for p in G.vs:
# 	if p['type'] == 'P':
# 		print G.neighbors(p)
# 		print G.incident(p)
# 		for n in G.neighbors(p):
# 			print G.es[ G.get_eid(p.index,n) ]['J'] 
# 			if not G.es[G.get_eid(p.index,n)]['J']:	
# 				print 'Dupa'
# 				G.es[G.get_eid(p.index,n)]['J'] = jMax
# 				jMax += 1
# 			print
# print G.es['J']

# for G in subproblems:
# 	print G






# layout = BFG.layout("kk") 				# Makes a full graph
# BFG.vs['label'] = BFG.vs['type']
# color_dict = {'F': 'green', 'P': 'blue', 'E': 'red'}
# BFG.vs['color'] = [ color_dict[ge] for ge in BFG.vs['type']]
# ig.plot(BFG, layout = layout, target=savePath+'/fullG.pdf')#,  bbox = (300, 300), margin = 20)


# ig.plot(BFG, layout = layout, target=savePath+'/suppByExp.pdf') 			# Severed bonds with experimentally unsupported nodes.

# color_dict = {'F': 'green', 'P': 'blue', 'E': 'red'}
# for i,G in enumerate(subproblems):
# 	layout = G.layout("kk")	
# 	G.vs['label'] = G.vs['type']
# 	G.vs['color'] = [ color_dict[ge] for ge in G.vs['type']]
# 	ig.plot( G, layout = layout, target= savePath + '/graph' + str(i) + '.pdf' )
