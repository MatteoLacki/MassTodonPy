##### The (ever growing) list of dependencies
from 	data.proteins import deconvolutionProblem, substanceP
import 	intervaltree as it
# from exampleData import theory, tol, minExpSupp, empiria, savePath
from 	collections import Counter, defaultdict
import 	igraph as ig
from 	scipy.optimize import linprog
import 	numpy as np
from 	math import sqrt

def genGraph(EXP, TH, tol, minExpSupp, realValuesTmp):
	tolIntervals = it.IntervalTree()
	V = []
	E = []
	for famNo in xrange(len(TH)):
		familyName = 'F_'+str(famNo)
		V.append({ 'name':familyName, 'potentialSupport': 0.0, 'type':'F' })
		for peakNo, info in enumerate(TH[famNo]):
			m, p = info
			peakName = 'P_'+str(peakNo)+'_'+str(famNo)
			V.append({ 'name':peakName, 'mass':m, 'intensity':p, 'type':'P' })
			E.append({ 'source': familyName, 'target': peakName})
			tolIntervals[ m-tol : m+tol ] = peakName
	realValues = {}
	for famNo in xrange(len(TH)):
		realValues['F_'+str(famNo)]=realValuesTmp[famNo]
	gNo = 0
	G 	= defaultdict( lambda: None )
	Gd 	= {}
	atheoretic = defaultdict(list)
	for eNo, info in enumerate(EXP):
		m, i = info
		intervals = tolIntervals[m]
		if intervals:
			eName = 'E_'+str(eNo)
			V.append({ 'name':eName, 'mass':m, 'intensity':i, 'type':'E' })
			grandpas = frozenset( interval.data for interval in intervals)
			if not G[grandpas]: 			# New G node.
				gName 	= 'G_'+str(gNo)
				gNo 	+= 1
				g = {'name':gName, 'intensity':0.0, 'type':'G'}
				V.append(g)
				Gd[gName] = g
				for grandpa in grandpas:
					E.append({'source':grandpa, 'target':gName})
				G[grandpas] = gName
			gName = G[grandpas]
			E.append({'source':eName, 'target':gName})
			Gd[gName]['intensity'] += i
		else:
			atheoretic['experimentalPeaks'].append(info)
	BFG = ig.Graph.DictList(vertices=V, edges=E)
	# This can be done using the BFS iterator. Simplify later..
	for f in BFG.vs(type='F'): 								# Establishing experimental support for individual compounds.
		for peakNo in BFG.neighbors(f): 					# Just the envelopes here. Root not yet added
			if BFG.neighborhood_size( peakNo, order=1 ) > 2: 	# Ascertain that there are some experimental groups G around the P peak. Not include P and F (order at most 1)
				f['potentialSupport'] += BFG.vs[peakNo]['intensity']
		if f['potentialSupport'] < minExpSupp:
			for peakNo in BFG.neighbors(f):
		 		BFG.delete_edges( 	BFG.get_eid(peakNo, nodeNo) for nodeNo in BFG.neighbors( peakNo )
		 							if nodeNo != f.index ) # Severe edges between peaks and experimental peaks, not between peaks and family nodes
	nodesToDelete = [] 						# Deletion of the G nodes without any envelope peaks around left
	for g in BFG.vs(type='G'):
		notNeighboringP = True	 			# Guard for absence of P
		neighborsNos = BFG.neighbors(g)
		while notNeighboringP and len(neighborsNos)>0:
			neighborNo 	= neighborsNos.pop()
			neighbor 	= BFG.vs[ neighborNo ]
			if neighbor['type'] == 'P':
				notNeighboringP = False
		if notNeighboringP:
			nodesToDelete.append(g.index)
	BFG.delete_vertices( nodesToDelete )
	subproblems = [ subG for subG in BFG.decompose() if {'F','E'} <= set( subG.vs['type']) ]
	return (realValues, subproblems)

def notRoot(tup, RootNo):
	a, b = tup
	if a == RootNo:
		return b
	else:
		return a

def anal(subG):
	subG.es['varNo'] = None
	subG.add_vertex( type='R',name='Root') 	# Add the root: its edges store the mixture variables (alphas) in the optimisation problem
	RootNo = subG.vs.find(name='Root').index
	# Optimise
	for varNo, f in enumerate( subG.vs( type_eq='F' ) ):
		subG.add_edge( 'Root', f['name'], varNo=varNo, varType='alpha' ) 	# Number the alphas
	# plott(subG, bbox = (1500, 1500), target='/Volumes/doom/Users/matteo/Dropbox/Science/MassSpectrometry/MassTodon/Visual/subG.pdf')
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
	estimates = [ (subG.vs[notRoot(e.tuple, RootNo)]['name'], res.x[e['varNo']]) for e in subG.es(varType='alpha') ]
	return estimates


def analyze(massNoise=0.01):
	tol 		= .05
	minExpSupp 	= .6
	chosenValues = [[10000, 20000, 15000],[5000, 23000, 6000]]
	realValuesTmp = [ v for cnt in chosenValues for v in cnt]
	EXP, TH 	= deconvolutionProblem( substanceP, chosenValues[0], massNoise=massNoise, digit=2)
	EXP2, TH2 	= deconvolutionProblem( substanceP+'A', chosenValues[1], massNoise=massNoise, digit=2)
	EXP.extend(EXP2)
	TH.extend(TH2)
	realValues, subproblems = genGraph(EXP, TH, tol, minExpSupp, realValuesTmp)
	results = [ Counter(dict(anal(G))) for G in subproblems ]
	I = iter(results)
	R = I.next()
	for i in I:
		R += i
	TotalError 		= sum( abs(R[fNo]-realValues[fNo]) for fNo in realValues )
	normalisation 	= sum( realValues[fNo] for fNo in realValues )
	relativeError 	= TotalError/normalisation
	return (TotalError, relativeError, normalisation)

Res = [ analyze(massNoise=noise) for noise in np.linspace(.01, .05, num=100) ]
RelErrors = [relativeError for TotalError, relativeError, normalisation in Res]


import pandas as pd
Results = pd.DataFrame( (noise, err) for noise, err in zip(np.linspace(.01, .05, num=100), RelErrors) )
Results.columns = ['peak_position_st_dev', 'error']

Results.to_csv(path_or_buf='/Volumes/doom/Users/matteo/Dropbox/Science/MassSpectrometry/MassTodon/Posters/R/anal.csv',index=False)

print Results
# def plott( G, **kwds ):
# 	color 	= {'F':'green', 'E':'pink', 'G':'red', 'P':'blue', 'R': 'yellow'}
# 	visual_style= {}
# 	visual_style['vertex_color']= [ color[c] for c in G.vs['type'] ]
# 	visual_style['vertex_label']= G.vs['name']
# 	visual_style.update(kwds)
# 	ig.plot( G, **visual_style )

# layout = BFG.layout_lgl()
# plott(BFG, layout=layout, bbox = (1500, 1500))
# plott(BFG, bbox = (1500, 1500), target='/Volumes/doom/Users/matteo/Dropbox/Science/MassSpectrometry/MassTodon/Visual/insilicoGraph.pdf')
