import intervaltree as it
from exampleData import theory, tol, minExpSupp, empiria 
from collections import Counter, defaultdict
import igraph as ig

BFG = ig.Graph(directed=True) # the BIG FUCKING GRAPH
tolIntervals = it.IntervalTree()

for famNo, spectrum in enumerate(theory):
	famName = 'F_' + str(famNo)
	BFG.add_vertex( type='F', name=famName, experimentalSupport=0.0 )
	for peakNo, info in enumerate(spectrum):
		mass, intensity = info
		peakName = 'P_' + str(peakNo) + '_' + str(famNo)
		BFG.add_vertex( type='P',name=peakName, mass=mass, intensity=intensity )
		BFG.add_edge( famName, peakName )
		tolIntervals[ mass-tol : mass+tol ] = peakName

for expPeakNo, info in enumerate(empiria):
	mass, intensity = info
	expPeakName = 'E_'+str(expPeakNo)
	BFG.add_vertex( type='E', name=expPeakName )
	for interval in tolIntervals[mass]:
		peakName = interval.data
		BFG.add_edge( expPeakName, peakName )
		BFG.add_edge( peakName, expPeakName )

	# Establishing experimental support for individual compounds.	
for family in BFG.vs(type='F'):
	for peakNo in BFG.neighborhood( family, order=1, mode='out' ):
		if peakNo != family.index:
				# BELOW: > 1 for the same peak is always included, 
				# as order is set at most 1, and we want exactly 1.1
			if BFG.neighborhood_size( peakNo, order=1, mode='out' ) > 1:
				family['experimentalSupport'] += BFG.vs[peakNo]['intensity']

poorSupport = [ BFG.neighborhood( family, order=1, mode='out' ) for family in BFG.vs( type='F', experimentalSupport_lt=minExpSupp )]
poorSupport = [ i for f in poorSupport for i in f ]

	## To make full graph
# layout = BFG.layout("kk")
# BFG.vs['label'] = BFG.vs['type']
# color_dict = {'F': 'green', 'P': 'blue', 'E': 'red'}
# BFG.vs['color'] = [ color_dict[ge] for ge in BFG.vs['type']]
# ig.plot(BFG, layout = layout, target='/Users/matteo/fullGraph.pdf')#,  bbox = (300, 300), margin = 20)

BFG.delete_vertices( poorSupport )	
	## After cleaving the experimentally unsupported part.
# layout = BFG.layout("kk")
# ig.plot(BFG, layout = layout, target='/Users/matteo/experimentallySupported.pdf')#,  bbox = (300, 300), margin = 20)

print BFG.clusters()
