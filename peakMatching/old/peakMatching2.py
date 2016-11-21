import intervaltree as iTree
from exampleData import theory, tol, totProbThr, empiria 
from collections import Counter, defaultdict
from itertools import combinations

tree = iTree.IntervalTree()
for famNo, spectrum in enumerate(theory):
	for peakNo, peak in enumerate(spectrum):
		mass, intensity = peak
		tree[ mass-tol : mass+tol ] = (peakNo, famNo)

	# removing families based on the inability to explain the empirical spectrum
notVisited 	= defaultdict(lambda : True)
totalProb 	= Counter()
for mass, intensity in empiria:
	for interval in tree[mass]:
		peakNo, famNo = interval.data	
		theoryMass, theoryIntensity = theory[famNo][peakNo]
		if notVisited[ (peakNo,famNo) ]:
			notVisited[ (peakNo,famNo) ] = False
			totalProb[famNo] += theoryIntensity

	# deleting intervals corresponding to the unwanted families
for interval in tree.copy():
	peakNo, famNo = interval.data
	if totalProb[famNo] < totProbThr:
		tree.remove(interval)

	# attributing experimental peaks to base intervals (that make up sums of tolerance intervals)
atheoretic = []
baseIntervals = Counter()

V = set()				# The verteces of the family numbers graph
E = defaultdict(set) 	# Adjecency table of the same graph

for mass, intensity in empiria:
	intervals = tree[mass]
	if not intervals:
		atheoretic.append( (mass,intensity) )
	else:
		familyClique = set([ interval.data[1] for interval in intervals ])
			# Add edges to the graph	
		for A, B in combinations( familyClique, 2 ):
			E[A].add(B)
			E[B].add(A)
			# Add vertices to the graph
		V |= familyClique
			# Builds up the intensity list.
		baseIntervals[ tuple([ interval.data for interval in intervals ]) ] += intensity	

# print baseIntervals 

	# Elucidating what defines the deconvolution subproblems
toVisit = set()
CCs = []

while V:
	if not toVisit:
		famNo = V.pop()
		CCs.append( set([famNo]) )
		CC = CCs[-1]
	else:
		famNo = toVisit.pop()
	for neighbour in E[famNo]: 
		if neighbour in V and not neighbour in toVisit:
			CC.add(neighbour)
			V.discard(neighbour)
			toVisit.add(neighbour)

# print baseIntervals.keys()
print CCs

cluster = []
clusters = [cluster]
for experimentPeaks in baseIntervals.keys():
	# print experimentPeaks
	for peakNo, famNo in experimentPeaks:
		cluster


print clusters
# for familiesInCluster in CCs:
# 	while familiesInCluster:
# 		famNo = familiesInCluster.pop()
# 		print famNo	

# 	clusters.append(cluster)
	# pass
	# clusters.append(  )




# for indexSet in baseIntervals.keys():
# 	__, famNo = indexSet[0]
