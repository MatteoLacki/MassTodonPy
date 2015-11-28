##### The (ever growing) list of dependencies
import intervaltree as it
from exampleData import theory, tol, minExpSupp, empiria, savePath
from collections import Counter, defaultdict
import igraph as ig
import pprint as PP
pp = PP.PrettyPrinter()
from scipy.optimize import linprog

##### The programme
BFG = ig.Graph() 											# the BIG FUCKING GRAPH (le grand graphe de merde)
tolIntervals = it.IntervalTree()

for famNo, spectrum in enumerate(theory):					# Insert family (F) and envelope (P) peaks in the BFG
	famName = 'F_' + str(famNo)
	BFG.add_vertex( type='F', name=famName, potentialSupport=0.0 )
	for peakNo, info in enumerate(spectrum):
		mass, intensity = info
		peakName = 'P_' + str(peakNo) + '_' + str(famNo)
		BFG.add_vertex( type='P',name=peakName, mass=mass, intensity=intensity )
		BFG.add_edge( famName, peakName )					# Add edges between P peaks belonging to an F node
		tolIntervals[ mass-tol : mass+tol ] = peakName


EGNo = 0
EG = defaultdict( lambda: None )
for eNo, info in enumerate(empiria):
	mass, intensity = info
	intervals = tolIntervals[mass]

	if intervals:
		Ename = 'E_' + str(eNo)
		BFG.add_vertex( type='E', name=Ename, mass=mass, intensity=intensity )
		grandparentNames = frozenset( interval.data for interval in intervals)
		
		if EG[ grandparentNames ] == None: 									# grouping node (G) with grandparents not yet existing
			EGname = 'G_' + str(EGNo)
			EGNo += 1
			BFG.add_vertex( type='G', name=EGname, intensity=0.0 ) 			# creates parent with appropriate grandparents
			for grandpaName in grandparentNames:
				BFG.add_edge( EGname, grandpaName ) 						# links parent with grandparents 
			EG[ grandparentNames ] = EGname 								# stores the parent name based on its grandparents
		
		EGname = EG[ grandparentNames ]
		BFG.add_edge( Ename, EGname )
		BFG.vs.find(name=EGname)['intensity'] += intensity
	else:
		BFG.add_vertex( type='E', name='atheoretic', mass=mass, intensity=intensity )	# a peak not within the theoretically predicted chemical compounds

def getAtheoreticalPeaks(G):
	'Masses and intensities of peaks that cannot be attributed to any theoretically observable chemical compound.'
	return [ 	( experim['mass'], experim['intensity'] ) for 
				experim in BFG.vs(type='E') if experim.degree() == 0 	]	

def getUnsupportedFamilies(G, minimalProbability):
	'Info on the chemical compounds with potentially less support in data than the minimalProbability value.'
	print 'IMPLEMENT IT'
	pass

for family in BFG.vs(type='F'): 							# Establishing experimental support for individual compounds.	
	for peakNo in BFG.neighbors(family): 					# Just the envelopes here. Root not yet added
		if BFG.neighborhood_size( peakNo, order=1 ) > 2: 	# Ascertain that there are some experimental groups G around the P peak. Not include P and F (order at most 1)
			family['potentialSupport'] += BFG.vs[peakNo]['intensity'] 	
	if family['potentialSupport'] < minExpSupp:
		for peakNo in BFG.neighbors(family):
	 		BFG.delete_edges( 	BFG.get_eid(peakNo, nodeNo) for nodeNo in BFG.neighbors( peakNo )  
	 							if nodeNo != family.index ) # Severe edges between peaks and experimental peaks, not between peaks and family nodes

nodesToDelete = [] 						# Deletion of the G nodes without any envelope peaks around left
for G in BFG.vs(type='G'):
	notNeighboringP = True	 			# Guard for absence of P
	neighborsNos = BFG.neighbors(G)
	while notNeighboringP and len(neighborsNos)>0:
		neighborNo 	= neighborsNos.pop()
		neighbor 	= BFG.vs[ neighborNo ]
		if neighbor['type'] == 'P':
			notNeighboringP = False
	if notNeighboringP:
		nodesToDelete.append(G.index)
BFG.delete_vertices( nodesToDelete )

subproblems = [ subG for subG in BFG.decompose() 				# subproblems for the deconvolution
				if {'F','E'} <= set( subG.vs['type']) ]		# Only good subproblems have both a family node and an empirical node 

for s in subproblems:
	print s

subG = subproblems[3]
subG.es['varNo'] = None

subG.add_vertex( type='R',name='Root') 	# Add the root: its edges store the mixture variables (alphas) in the optimisation problem
RootNo = subG.vs.find(name='Root').index
for varNo, f in enumerate( subG.vs( type_eq='F' ) ):
	subG.add_edge( 'Root', f['name'], varNo=varNo, varType='alpha' ) 	# Number the alphas
varNoMax = varNo+1

AeqSparse = []	# A scarse matrix representation of the equality constraints in the optimisation problem
for i,P in enumerate( subG.vs( type = 'P' ) ):
	if len( subG.neighbors(P) ) > 1 :			# If no G nodes around, should not create an equality constraint
		for nodeNo in subG.neighbors(P):
			node = subG.vs[nodeNo]
			isFamilyNode = node['type'] == 'F'
			if isFamilyNode: 								# edge between either:
				edge = subG.es[ subG.get_eid(nodeNo, RootNo) ] 		# the root and a family 
			else:													# or
				edge = subG.es[ subG.get_eid( P.index, nodeNo) ]		# an experimental peak and an envelope peak
			
			if edge['varNo'] == None: 			
				edge['varNo'] 	= varNoMax 	# Assign number to the flow variable in the optimisation problem
				edge['varType'] = 'flow' 	# Tag as flow
				varNoMax += 1
			
			if isFamilyNode:
				AeqSparse.append(( i, edge['varNo'], -float(P['intensity']) ))
			else:
				AeqSparse.append(( i, edge['varNo'], 1.0 ))

normalize = False
if normalize:
	for alphaEdge in subG.es[ subG.adjacent('Root') ]: 		# This assures that alphas sum to one. This is not necessary or wanted...
		AeqSparse.append((i,alphaEdge['varNo'],1.0))
	
zeros = [0.0] * varNoMax 	# Densify the sparse representation
Aeq = []
I = i+1
for i in xrange(I):
	Aeq.append(zeros[:])
for i,j,v in AeqSparse:
	Aeq[i][j] = v
Beq = [0.0] * len(Aeq)

AineqSparse = [] 	# Prepare a scarse matrix representation of the inequality constraints in the optimisation problem
Bineq = []			# Not really sparse here. But can make it later...

for i, G in enumerate( subG.vs( type = 'G' ) ):
	Bineq.append( G['intensity'] )
	for nodeNo in subG.neighbors(G):
		node = subG.vs[nodeNo]
		if node['type'] == 'P':
			e = subG.es[ subG.get_eid( nodeNo, G.index ) ]
			AineqSparse.append(( i, e['varNo'], 1.0 ))

I = i+1
Aineq = []
for i in xrange(I):
	Aineq.append(zeros[:])
for i,j,v in AineqSparse:
	Aineq[i][j] = v 

c0Sparse = [ (e['varNo'], 1.) for e in subG.es( varType = 'flow' ) ]
c0 = [0.0]*varNoMax 
for i,v in c0Sparse:
	c0[i]=-v 		# minus for the optimisation is minimising and we want max flow.

c1 = c0[:]

res = linprog(c0, A_ub=Aineq, b_ub=Bineq, options={'disp':True}, A_eq=Aeq, b_eq=Beq )
print res.x

# # pp.pprint( Aineq )
# # print Bineq

# # print [ a for a in G.es(varNo_ne=None )]
# # print [ a for a in G.es['varNo']  ]

# # print [ (a['mass'],a['intensity']) for a in G.vs(type='E')]


layout = subG.layout("kk") 				# Makes a full graph
subG.vs['label'] = subG.vs['type']
color_dict = {'F': 'green', 'P': 'blue', 'E': 'pink', 'G': 'red', 'R': 'yellow' }
subG.vs['color'] = [ color_dict[ge] for ge in subG.vs['type']]
subG.es['label'] = subG.es['varNo']
ig.plot(subG, layout = layout, target=savePath+'/artificialNode.pdf')#,  bbox = (300, 300), margin = 20)


# layout = BFG.layout("kk") 				# Makes a full graph
# BFG.vs['label'] = BFG.vs['type']
# color_dict = {'F': 'green', 'P': 'blue', 'E': 'pink', 'G': 'red', 'R': 'yellow' }
# BFG.vs['color'] = [ color_dict[ge] for ge in BFG.vs['type']]
# ig.plot(BFG, layout = layout, target=savePath+'/suppByExp.pdf')#,  bbox = (300, 300), margin = 20)



# # layout = G.layout("kk") 				# Makes a full graph
# # G.vs['label'] = G.vs['type']
# # color_dict = {'F': 'green', 'P': 'blue', 'E': 'red', 'R': 'yellow' }
# # G.vs['color'] = [ color_dict[ge] for ge in G.vs['type']]
# # G.es['label'] = G.es['varNo']
# # ig.plot(G, layout = layout, target=savePath+'/artificialNode.pdf')#,  bbox = (300, 300), margin = 20)


# # layout = BFG.layout("kk") 				# Makes a full graph
# # BFG.vs['label'] = BFG.vs['type']
# # color_dict = {'F': 'green', 'P': 'blue', 'E': 'red'}
# # BFG.vs['color'] = [ color_dict[ge] for ge in BFG.vs['type']]
# # ig.plot(BFG, layout = layout, target=savePath+'/fullG2.pdf')#,  bbox = (300, 300), margin = 20)


# # # ig.plot(BFG, layout = layout, target=savePath+'/suppByExp.pdf') 			# Severed bonds with experimentally unsupported nodes.

# # # color_dict = {'F': 'green', 'P': 'blue', 'E': 'red'}
# # # for i,G in enumerate(subproblems):
# # # 	layout = G.layout("kk")	
# # # 	G.vs['label'] = G.vs['type']
# # # 	G.vs['color'] = [ color_dict[ge] for ge in G.vs['type']]
# # # 	ig.plot( G, layout = layout, target= savePath + '/graph' + str(i) + '.pdf' )
