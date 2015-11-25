from collections import defaultdict

V = set([1,2,3,4,5,6,7,8,9,10])
E = defaultdict(set) 

	# Sets are needed because of the construction of the graph. 
	# Not indispensable, but otherwise we would need ... a table with translations 
	# or just one renaming.
E[1] = set([5])
E[2] = set([5])
E[3] = set([5])
E[4] = set([7])
E[5] = set([1,2,3])
E[6] = set([7])
E[7] = set([4,6])
E[8] = set([9,10])
E[9] = set([8,10])
E[10] = set([8,9])

# This is in fact something between BFS and DFS: the queue is not FIFO nor LIFO
# But it works.
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
			# neighbour in V - not yet visited
			# not neighbour in toVisit - not on the toVisit list
		if neighbour in V and not neighbour in toVisit:
			CC.add(neighbour)
			V.discard(neighbour)
			toVisit.add(neighbour)
	
print CCs
