import intervaltree as it
from exampleData import theory, tol, minExpSupp, empiria, savePath
from collections import Counter, defaultdict
import igraph as ig

BFG = ig.Graph() 							# the BIG FUCKING GRAPH
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
	pass

for family in BFG.vs(type='F'): 										# Establishing experimental support for individual compounds.	
	for peakNo in BFG.neighbors(family):
		if BFG.neighborhood_size( peakNo, order=1 ) > 2: 				# > 2 for the PEAK and FAMILY peak are within NEIGBORHOOD( ORDER = 1 )
			family['potentialSupport'] += BFG.vs[peakNo]['intensity'] 	
	if family['potentialSupport'] < minExpSupp:
		for peakNo in BFG.neighbors(family):
	 		BFG.delete_edges( 
	 			BFG.get_eid(peakNo, nodeNo) for nodeNo in BFG.neighbors( peakNo ) 
	 			if nodeNo != family.index ) 			# Severe edges between peaks and experimental peaks, but not between peaks and family node

subproblems = [ G for G in BFG.decompose() 
				if {'F','E'} <= set(G.vs['type']) ]		# Only good subproblems have both a family node and an empirical node 

G = subproblems[1]
print G


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
